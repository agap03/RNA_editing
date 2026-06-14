#!/usr/bin/env python3
"""
Survival analysis for a single cancer type, stratified by high vs low AEI.

Merges one TCGA survival TSV with tumor AEI from REDIportal, splits patients
at the median AEI within that cancer type, and runs Kaplan-Meier analysis
with a log-rank test.
"""
import argparse
import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from plot_utils import save_figure, significance_stars
from redi_utils import derive_tumor_cancer_type

SCRIPT = "survival_analysis"

# Minimum cohort sizes required to run a stable log-rank comparison.
MIN_GROUP_SIZE = 5
MIN_TOTAL = 10


def tcga_patient_id(sample):
    """Reduce a full TCGA aliquot barcode to the 3-part patient identifier."""
    parts = str(sample).split("-")
    if len(parts) >= 3:
        return "-".join(parts[:3])
    return str(sample)


def cancer_type_from_clinical_path(path):
    """Parse BLCA from a file named TCGA-BLCA.survival.tsv."""
    match = re.search(r"TCGA-([A-Z0-9]+)\.survival", os.path.basename(path))
    if match is None:
        raise ValueError(f"Could not parse cancer type from clinical file: {path}")
    return match.group(1)


def load_clinical_file(clinical_path, cancer_type):
    """Load one TCGA survival table with OS time (months) and event indicator."""
    clinical = pd.read_csv(clinical_path, sep="\t")
    clinical["Cancer_type"] = cancer_type
    clinical["OS.time"] = pd.to_numeric(clinical["OS.time"], errors="coerce")
    clinical["OS"] = pd.to_numeric(clinical["OS"], errors="coerce")
    return clinical.dropna(subset=["OS.time", "OS", "_PATIENT"])


def load_tumor_aei(cancers_path, cancer_type):
    """
    Load tumor-only AEI for one cancer type.

    Multiple aliquots from the same patient are collapsed to the median AEI.
    """
    tumor = pd.read_csv(cancers_path)
    tumor = derive_tumor_cancer_type(tumor)
    tumor = tumor[(tumor["Status"] == "tumor") & (tumor["Cancer_type"] == cancer_type)].copy()
    tumor["_PATIENT"] = tumor["Sample"].apply(tcga_patient_id)
    tumor["AEI"] = pd.to_numeric(tumor["AEI"], errors="coerce")
    tumor = tumor.dropna(subset=["AEI", "_PATIENT"])
    tumor = tumor[(tumor["AEI"].isna()) | (tumor["AEI"] >= 0)]
    return tumor.groupby("_PATIENT", as_index=False)["AEI"].median()


def assign_aei_groups(df):
    """Split patients into high vs low AEI using the within-cancer median."""
    median_aei = df["AEI"].median()
    df = df.copy()
    df["AEI_group"] = df["AEI"].apply(
        lambda value: "High AEI" if value >= median_aei else "Low AEI"
    )
    df["AEI_median"] = median_aei
    return df


def median_survival(df, group):
    """Median survival time for one AEI group, if estimable."""
    subset = df[df["AEI_group"] == group]
    kmf = KaplanMeierFitter()
    kmf.fit(subset["OS.time"], subset["OS"])
    return round(kmf.median_survival_time_, 1) if kmf.median_survival_time_ is not None else None


def run_survival_analysis(df, cancer_type, output_dir):
    """Fit Kaplan-Meier curves and run a two-group log-rank test."""
    low = df[df["AEI_group"] == "Low AEI"]
    high = df[df["AEI_group"] == "High AEI"]

    if len(df) < MIN_TOTAL or len(low) < MIN_GROUP_SIZE or len(high) < MIN_GROUP_SIZE:
        print(f"[{SCRIPT}] Skipped: insufficient samples for {cancer_type}")
        return None

    result = logrank_test(
        low["OS.time"], high["OS.time"],
        low["OS"], high["OS"],
    )

    fig, ax = plt.subplots(figsize=(7, 5))
    kmf = KaplanMeierFitter()
    for group, color in [("Low AEI", "steelblue"), ("High AEI", "salmon")]:
        subset = df[df["AEI_group"] == group]
        kmf.fit(subset["OS.time"], subset["OS"], label=f"{group} (n={len(subset)})")
        kmf.plot(ax=ax, ci_show=True, color=color)

    sig = significance_stars(result.p_value)
    ax.set_title(
        f"{cancer_type}: overall survival by AEI group\n"
        f"log-rank p = {result.p_value:.4g} {sig}",
        fontsize=12,
    )
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Overall survival probability")
    ax.set_ylim(0, 1.05)
    plt.tight_layout()
    save_figure(SCRIPT, os.path.join(output_dir, f"{cancer_type}_km.png"))

    print(
        f"[{SCRIPT}] {cancer_type}: "
        f"low n={len(low)}, high n={len(high)}, p={result.p_value:.4g} {sig}"
    )

    return {
        "Cancer_type": cancer_type,
        "n_total": len(df),
        "n_low_AEI": len(low),
        "n_high_AEI": len(high),
        "events_low_AEI": int(low["OS"].sum()),
        "events_high_AEI": int(high["OS"].sum()),
        "median_AEI_cutoff": round(df["AEI_median"].iloc[0], 3),
        "median_survival_low_AEI": median_survival(df, "Low AEI"),
        "median_survival_high_AEI": median_survival(df, "High AEI"),
        "logrank_p_value": round(result.p_value, 6),
        "Significance": sig,
    }


def main():
    parser = argparse.ArgumentParser(description="Survival analysis by high vs low AEI for one cancer type")
    parser.add_argument("--cancers", "-c", required=True, help="Path to cancer_REDIportal.csv")
    parser.add_argument("--clinical", required=True, help="Path to one TCGA survival TSV file")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    args = parser.parse_args()

    cancer_type = cancer_type_from_clinical_path(args.clinical)
    os.makedirs(args.output, exist_ok=True)

    print(f"[{SCRIPT}] Cancer type: {cancer_type}")
    print(f"[{SCRIPT}] Loading clinical survival data...")
    clinical = load_clinical_file(args.clinical, cancer_type)

    print(f"[{SCRIPT}] Loading tumor AEI data...")
    aei = load_tumor_aei(args.cancers, cancer_type)

    # Match clinical and molecular data at the TCGA patient level.
    merged = clinical.merge(aei, on="_PATIENT", how="inner")
    merged = merged.dropna(subset=["OS.time", "OS", "AEI"])
    merged = assign_aei_groups(merged)

    merged_path = os.path.join(args.output, f"{cancer_type}_merged.csv")
    merged.to_csv(merged_path, index=False)
    print(f"[{SCRIPT}] Saved: {os.path.abspath(merged_path)}")

    result = run_survival_analysis(merged, cancer_type, args.output)
    logrank_path = os.path.join(args.output, f"{cancer_type}_logrank.csv")
    km_path = os.path.join(args.output, f"{cancer_type}_km.png")

    # Snakemake expects all declared outputs; write placeholders when skipped.
    if result is None:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            f"{cancer_type}: insufficient samples for survival analysis",
            ha="center", va="center", fontsize=12,
        )
        plt.tight_layout()
        save_figure(SCRIPT, km_path)
        pd.DataFrame([{
            "Cancer_type": cancer_type,
            "Status": "insufficient_samples",
        }]).to_csv(logrank_path, index=False)
    else:
        pd.DataFrame([result]).to_csv(logrank_path, index=False)

    print(f"[{SCRIPT}] Saved: {os.path.abspath(logrank_path)}")
    print(f"[{SCRIPT}] Done. Output directory: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
