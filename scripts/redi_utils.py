"""Shared utilities for loading and preparing REDIportal data."""

import re

import pandas as pd

# Maps TCGA cancer codes to keywords found in normal-tissue sample names.
NORMAL_TISSUE_MAP = {
    "ACC": ["Adrenal Gland"],
    "BLCA": ["Bladder"],
    "BRCA": ["Breast - Mammary"],
    "CESC": ["Cervix -"],
    "CHOL": ["Liver"],
    "COAD": ["Colon -", "Small Intestine -"],
    "READ": ["Colon -", "Small Intestine -"],
    "ESCA": ["Esophagus -"],
    "GBM": ["Brain -"],
    "LGG": ["Brain -"],
    "HNSC": ["Salivary Gland"],
    "KICH": ["Kidney -"],
    "KIRC": ["Kidney -"],
    "KIRP": ["Kidney -"],
    "LIHC": ["Liver"],
    "LUAD": ["Lung"],
    "LUSC": ["Lung"],
    "OV": ["Ovary"],
    "PAAD": ["Pancreas"],
    "PCPG": ["Adrenal Gland"],
    "PRAD": ["Prostate"],
    "SKCM": ["Skin -"],
    "STAD": ["Stomach"],
    "TGCT": ["Testis"],
    "THCA": ["Thyroid"],
    "UCEC": ["Uterus"],
    "UCS": ["Uterus"],
}


def tissue_matches_keyword(tissue, keyword):
    """Match a GTEx tissue label to a map keyword without substring false positives."""
    tissue_str = str(tissue)
    if keyword.startswith("^") or any(character in keyword for character in "[]()$+"):
        return bool(re.search(keyword, tissue_str, flags=re.IGNORECASE))
    if " -" not in keyword and " " not in keyword:
        return tissue_str.lower() == keyword.lower()
    return keyword.lower() in tissue_str.lower()


def tissue_matches_keywords(tissue, keywords):
    """Return True when any keyword matches the tissue label."""
    return any(tissue_matches_keyword(tissue, keyword) for keyword in keywords)


def match_tissue_to_cancer_type(tissue):
    """Match a normal tissue name to a cancer type using NORMAL_TISSUE_MAP."""
    for cancer_type, keywords in NORMAL_TISSUE_MAP.items():
        if tissue_matches_keywords(tissue, keywords):
            return cancer_type
    return None


def derive_tumor_cancer_type(df):
    """Extract TCGA project code (e.g. BLCA) from tumor sample metadata."""
    df = df.copy()
    if "Body Site/Study" not in df.columns:
        return df
    # Primary pattern: TCGA-BLCA in the body site field.
    df["Cancer_type"] = (
        df["Body Site/Study"].astype(str)
        .str.extract(r"TCGA-([A-Za-z0-9_]+)", expand=False)
    )
    # Fallback for non-standard labels: use the last token after splitting on "-".
    df["Cancer_type"] = df["Cancer_type"].fillna(
        df["Body Site/Study"].astype(str).str.split("-").str[-1].str.strip()
    )
    return df


def derive_normal_cancer_type(df):
    """Assign a TCGA-like cancer code to normal samples via tissue keywords."""
    df = df.copy()
    if "Body Site/Study" in df.columns:
        df["Cancer_type"] = df["Body Site/Study"].apply(match_tissue_to_cancer_type)
    return df


def clean_labels(df, col="Cancer_type"):
    """Strip the TCGA- prefix for compact plot labels."""
    df = df.copy()
    df["Label"] = df[col].str.replace("TCGA-", "", regex=False)
    return df


def _coerce_numeric_columns(df):
    """Convert AEI and expression columns to numeric, coercing invalid values to NaN."""
    df = df.copy()
    for col in ["AEI", "ADAR", "ADARB1", "ADARB2"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _remove_negative_values(df, sample_kind, script=None):
    """Drop rows with negative AEI or expression values, which are invalid in REDIportal."""
    df = df.copy()
    for col in ["AEI", "ADAR", "ADARB1", "ADARB2"]:
        if col not in df.columns:
            continue
        before = len(df)
        df = df[(df[col].isna()) | (df[col] >= 0)]
        after = len(df)
        if script and before != after:
            print(f"[{script}] Removed {before - after} {sample_kind} rows with {col} < 0")
    return df


def load_redi_data(cancers_path, normal_path, script=None):
    """Load REDIportal CSVs with shared cancer-type derivation and cleaning."""
    tumor = derive_tumor_cancer_type(pd.read_csv(cancers_path))
    normal = derive_normal_cancer_type(pd.read_csv(normal_path))
    tumor = _coerce_numeric_columns(tumor)
    normal = _coerce_numeric_columns(normal)
    tumor = _remove_negative_values(tumor, "tumor", script)
    normal = _remove_negative_values(normal, "normal", script)
    return tumor, normal


def correct_aei_batch_effect(combined):
    """
    Correct systematic AEI offset between tumor (TCGA) and normal (GTEx) samples.

    Each batch is shifted so its median matches the global median of all matched
    tumor and normal AEI values combined.
    """
    combined = combined.copy()
    global_median = combined["AEI"].median()
    tumor_mask = combined["Status"] == "tumor"
    normal_mask = combined["Status"] == "normal"
    tumor_median = combined.loc[tumor_mask, "AEI"].median()
    normal_median = combined.loc[normal_mask, "AEI"].median()
    tumor_shift = global_median - tumor_median
    normal_shift = global_median - normal_median
    combined.loc[tumor_mask, "AEI"] = combined.loc[tumor_mask, "AEI"] + tumor_shift
    combined.loc[normal_mask, "AEI"] = combined.loc[normal_mask, "AEI"] + normal_shift
    return combined, {
        "global_median": global_median,
        "tumor_shift": tumor_shift,
        "normal_shift": normal_shift,
    }


def get_matched_tumor_normal_data(tumor, normal, normal_map=NORMAL_TISSUE_MAP, batch_correct=False):
    """
    Return tumor and normal rows only for cancer types with matched normal reference.

    Normal matching uses the same tissue-keyword logic as AEI_analysis.
    Optionally applies global median batch correction to tumor and normal AEI values.
    """
    tumor = tumor[tumor["Status"] == "tumor"].copy()
    normal = normal[normal["Status"] == "normal"].copy()

    tumor_parts = []
    normal_parts = []

    for cancer_type in sorted(tumor["Cancer_type"].dropna().unique()):
        keywords = normal_map.get(cancer_type, [])
        if not keywords:
            continue

        matched_normal = normal[
            normal["Body Site/Study"].apply(
                lambda tissue: tissue_matches_keywords(tissue, keywords)
            )
        ].copy()
        if matched_normal.empty:
            continue

        matched_normal["Cancer_type"] = cancer_type
        tumor_parts.append(tumor[tumor["Cancer_type"] == cancer_type])
        normal_parts.append(matched_normal)

    if not tumor_parts:
        return pd.DataFrame(), None

    combined = pd.concat(tumor_parts + normal_parts, ignore_index=True)
    batch_info = None
    if batch_correct:
        combined, batch_info = correct_aei_batch_effect(combined)
    return combined, batch_info
