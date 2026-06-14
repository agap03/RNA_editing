#!/usr/bin/env python3
"""
Integrate RNA-SBS1 signature activity with REDIportal AEI statistics.

Uses outputs from the REIA branch (signature activities) and the REDI branch
(AEI summaries and ADAR correlations) to produce cross-cohort scatter plots.
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from plot_utils import save_figure

SCRIPT = "combined_analysis_figures"


def main():
    parser = argparse.ArgumentParser(description="Analyze RNA-SBS1 activity relationships.")
    parser.add_argument("--activity", "-a", required=True, help="Path to Decompose_Solution_Activities.txt")
    parser.add_argument("--correlation", "-c", required=True, help="Path to AEI_vs_ADAR.csv")
    parser.add_argument("--aei", "-e", required=True, help="Path to AEI_summary_statistics.csv")
    parser.add_argument("--output", "-o", required=True, help="Output directory for figures")
    args = parser.parse_args()

    output_dir = args.output
    os.makedirs(output_dir, exist_ok=True)

    print(f"[{SCRIPT}] Loading activity, correlation, and AEI data...")

    activity = pd.read_csv(args.activity, sep="\t")
    correlation = pd.read_csv(args.correlation)
    aei = pd.read_csv(args.aei)

    # Column "0" is RNA-SBS1 activity in the SigProfiler activities file.
    activity["Cancer_type"] = activity["Samples"]
    activity = activity[["Cancer_type", "0"]].rename(columns={"0": "Activity"})

    # Figure 1: cancer-level RNA-SBS1 activity vs mean tumor AEI.
    df = pd.merge(aei, activity, on="Cancer_type", how="inner")

    plt.figure(figsize=(8, 6))
    plt.scatter(df["mean"], df["Activity"])
    for _, row in df.iterrows():
        plt.annotate(row["Cancer_type"], (row["mean"], row["Activity"]), fontsize=8)

    r, p = pearsonr(df["Activity"], df["mean"])
    rho, p_s = spearmanr(df["Activity"], df["mean"])
    plt.text(
        0.05, 0.95,
        f"Pearson r = {r:.3f}\nPearson p value = {p:.3g}\n"
        f"Spearman rho = {rho:.3f}\nSpearman p value = {p_s:.3g}",
        verticalalignment="top",
        transform=plt.gca().transAxes,
    )
    plt.xlabel("Mean AEI (%)")
    plt.ylabel("Activity of signature RNA-SBS1")
    plt.title("Activity vs Mean AEI")
    plt.tight_layout()
    save_figure(SCRIPT, os.path.join(output_dir, "activity_vs_mean_AEI.png"))

    # Figure 2: univariate view of RNA-SBS1 activity across cancer types.
    plt.figure(figsize=(8, 6))
    plt.bar(activity["Cancer_type"], activity["Activity"])
    plt.xticks(rotation=90)
    plt.xlabel("Cancer type")
    plt.ylabel("Activity")
    plt.title("Activity of signature RNA-SBS1 across cancer types")
    plt.tight_layout()
    save_figure(SCRIPT, os.path.join(output_dir, "activity_by_cancer_type.png"))

    # Figure 3: does higher signature activity track with stronger AEI-ADAR correlation?
    df = pd.merge(correlation, activity, on="Cancer_type", how="inner")

    plt.figure(figsize=(8, 6))
    plt.scatter(df["Activity"], df["Pearson_r"])
    for _, row in df.iterrows():
        plt.annotate(row["Cancer_type"], (row["Activity"], row["Pearson_r"]), fontsize=8)

    r, p = pearsonr(df["Activity"], df["Pearson_r"])
    rho, p_s = spearmanr(df["Activity"], df["Pearson_r"])
    plt.text(
        0.65, 0.95,
        f"Pearson r = {r:.3f}\nPearson p value = {p:.3g}\n"
        f"Spearman rho = {rho:.3f}\nSpearman p value = {p_s:.3g}",
        transform=plt.gca().transAxes,
        verticalalignment="top",
    )
    plt.xscale("log")
    plt.xlabel("Activity (log scale)")
    plt.ylabel("Pearson correlation (r)")
    plt.title("Activity vs Pearson correlation")
    plt.tight_layout()
    save_figure(SCRIPT, os.path.join(output_dir, "activity_vs_pearson_correlation.png"))

    print(f"[{SCRIPT}] Done. Output directory: {os.path.abspath(output_dir)}")


if __name__ == "__main__":
    main()
