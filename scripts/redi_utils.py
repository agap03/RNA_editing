"""Shared utilities for loading and preparing REDIportal data."""

import pandas as pd

# Maps TCGA cancer codes to keywords found in normal-tissue sample names.
NORMAL_TISSUE_MAP = {
    "BLCA": ["Bladder", "Urinary", "Urothelium"],
    "ESCA": ["Esophagus", "Esophageal"],
    "GBM": ["Brain", "Cereb", "Neural"],
    "LUAD": ["Lung", "Pulmonary", "Bronchus"],
    "LUSC": ["Lung", "Pulmonary", "Bronchus"],
    "BRCA": ["Breast", "Mammary"],
    "COAD": ["Colon", "Colorectal", "Intestine"],
    "HNSC": ["Head", "Neck", "Oral", "Throat"],
    "KIRC": ["Kidney", "Renal"],
    "LIHC": ["Liver", "Hepatic"],
    "OV": ["Ovary", "Ovarian"],
    "PAAD": ["Pancreas", "Pancreatic"],
    "PRAD": ["Prostate"],
    "SKCM": ["Skin", "Cutaneous", "Melanocyte"],
    "UCEC": ["Uterus", "Endometrium"],
}


def match_tissue_to_cancer_type(tissue):
    """Match a normal tissue name to a cancer type using NORMAL_TISSUE_MAP."""
    tissue_str = str(tissue).lower()
    for cancer_type, keywords in NORMAL_TISSUE_MAP.items():
        for keyword in keywords:
            if keyword.lower() in tissue_str:
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
