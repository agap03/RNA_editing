#!/usr/bin/env python3
"""
Generate REDIportal summary figures.

Visualises sample-level AEI distributions and expression correlations.
Statistical annotations for expression scatter plots come from AEI_analysis
outputs rather than being recomputed here.
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

from plot_utils import save_figure
from redi_utils import clean_labels, load_redi_data

SCRIPT = "REDI_statistics_figures"


def _correlation_stats_by_label(correlation_df):
    """Index precomputed per-cancer correlation stats by plot label."""
    if correlation_df is None or correlation_df.empty:
        return {}
    labeled = clean_labels(correlation_df)
    return labeled.set_index("Label").to_dict("index")


def figure1_aei_across_cancers(tumor, outdir):
    """Boxplot of AEI across all tumor cancer types."""
    tumor = clean_labels(tumor)
    order = sorted(tumor["Label"].unique())

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(data=tumor, x="Label", y="AEI", hue="Label",
                order=order, palette="Set2", ax=ax, legend=False,
                flierprops=dict(marker="o", markersize=2, alpha=0.3))
    sns.stripplot(data=tumor, x="Label", y="AEI",
                  order=order, color="black",
                  alpha=0.2, size=2, ax=ax, jitter=True)
    ax.set_title("ADAR editing activity (AEI) across cancer types", fontsize=13)
    ax.set_xlabel("Cancer type", fontsize=11)
    ax.set_ylabel("AEI (%)", fontsize=11)
    plt.xticks(rotation=90)
    plt.tight_layout()
    save_figure(SCRIPT, os.path.join(outdir, "AEI_across_cancers.png"))


def figure_aei_vs_expression(tumor, correlation_df, gene, outdir):
    """
    One scatter panel per cancer type for AEI vs gene expression.

    The regression line is drawn from sample-level data; the title statistics
    are read from the AEI_analysis CSV for consistency with the stats tables.
    """
    tumor = clean_labels(tumor)
    stats_by_label = _correlation_stats_by_label(correlation_df)
    labels = sorted(tumor["Label"].unique())
    n = len(labels)

    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4))
    if n == 1:
        axes = [axes]

    for ax, label in zip(axes, labels):
        subset = tumor[tumor["Label"] == label].dropna(subset=["AEI", gene])
        ax.scatter(subset[gene], subset["AEI"], alpha=0.4, s=8, color="steelblue")

        if len(subset) > 2:
            m, b, _, _, _ = stats.linregress(subset[gene], subset["AEI"])
            x_line = pd.Series([subset[gene].min(), subset[gene].max()])
            ax.plot(x_line, m * x_line + b, color="red", linewidth=1.5)

        if label in stats_by_label:
            row = stats_by_label[label]
            ax.set_title(
                f"{label}\nr={row['Pearson_r']:.2f} {row['Significance']}",
                fontsize=10,
            )
        else:
            ax.set_title(label, fontsize=10)

        ax.set_xlabel(f"{gene} (TPM)", fontsize=9)
        ax.set_ylabel("AEI (%)", fontsize=9)

    plt.suptitle(f"AEI vs {gene} expression per cancer type", fontsize=12, y=1.02)
    plt.tight_layout()
    save_figure(SCRIPT, os.path.join(outdir, f"AEI_vs_{gene}.png"), bbox_inches="tight")


def figure3_tumor_vs_normal(combined, outdir):
    """Side-by-side boxplots of AEI in tumor and matched normal samples."""
    combined = combined[combined["Cancer_type"].notna()].copy()
    combined = clean_labels(combined)
    order = sorted(combined["Label"].unique())

    fig, ax = plt.subplots(figsize=(15, 5))
    sns.boxplot(data=combined, x="Label", y="AEI",
                hue="Status", order=order,
                palette={"tumor": "salmon", "normal": "lightblue"},
                ax=ax,
                flierprops=dict(marker="o", markersize=2, alpha=0.3))
    plt.xticks(rotation=90)
    ax.set_title("AEI in tumor vs normal tissue", fontsize=13)
    ax.set_xlabel("Cancer type", fontsize=11)
    ax.set_ylabel("AEI (%)", fontsize=11)
    ax.legend(title="Status")
    plt.tight_layout()
    save_figure(SCRIPT, os.path.join(outdir, "tumor_vs_normal.png"))


def main():
    parser = argparse.ArgumentParser(description="Generate figures from REDIportal data")
    parser.add_argument("--cancers", "-c", required=True, help="Path to cancer_REDIportal.csv")
    parser.add_argument("--normal", "-n", required=True, help="Path to normal_REDIportal.csv")
    parser.add_argument("--correlations", required=True,
                        help="Directory with AEI_vs_<gene>.csv files from AEI_analysis")
    parser.add_argument("--output", "-o", required=True, help="Output directory for figures")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    print(f"[{SCRIPT}] Loading cancer and normal data...")
    tumor, normal = load_redi_data(args.cancers, args.normal, script=SCRIPT)

    combined = pd.concat([tumor, normal], ignore_index=True)
    tumor = combined[combined["Status"] == "tumor"].copy()

    if len(tumor) == 0:
        print(f"[{SCRIPT}] Skipped: no tumor samples found in input")
        return

    figure1_aei_across_cancers(tumor, args.output)

    for gene in ["ADAR", "ADARB1", "ADARB2"]:
        correlation_path = os.path.join(args.correlations, f"AEI_vs_{gene}.csv")
        if gene not in tumor.columns:
            print(f"[{SCRIPT}] Skipped: AEI vs {gene} figure ({gene} column missing)")
            continue
        if not os.path.exists(correlation_path):
            print(f"[{SCRIPT}] Skipped: AEI vs {gene} figure ({correlation_path} missing)")
            continue
        correlation_df = pd.read_csv(correlation_path)
        figure_aei_vs_expression(tumor, correlation_df, gene, args.output)

    figure3_tumor_vs_normal(combined, args.output)

    print(f"[{SCRIPT}] Done. Output directory: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
