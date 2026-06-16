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

from plot_utils import add_significance_legend, save_figure, significance_stars
from redi_utils import clean_labels, get_matched_tumor_normal_data, load_redi_data

SCRIPT = "REDI_statistics_figures"
EXPRESSION_GENES = ["ADAR", "ADARB1", "ADARB2"]


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


def _load_correlation_tables(correlations_dir):
    """Load per-gene AEI correlation tables produced by AEI_analysis."""
    tables = {}
    for gene in EXPRESSION_GENES:
        path = os.path.join(correlations_dir, f"AEI_vs_{gene}.csv")
        if os.path.exists(path):
            tables[gene] = pd.read_csv(path)
    return tables


def figure_aei_expression_heatmap(correlation_tables, outdir):
    """Heatmap of per-cancer Pearson r between AEI and ADAR family expression."""
    available_genes = [gene for gene in EXPRESSION_GENES if gene in correlation_tables]
    if not available_genes:
        print(f"[{SCRIPT}] Skipped: AEI expression correlation heatmap (no correlation tables)")
        return

    long_parts = []
    for gene in available_genes:
        part = correlation_tables[gene][["Cancer_type", "Pearson_r", "Significance"]].copy()
        part["Gene"] = gene
        long_parts.append(part)

    long_df = pd.concat(long_parts, ignore_index=True)
    matrix = long_df.pivot(index="Cancer_type", columns="Gene", values="Pearson_r")[available_genes]
    sig_matrix = long_df.pivot(index="Cancer_type", columns="Gene", values="Significance")[available_genes]

    if "ADAR" in matrix.columns:
        matrix = matrix.sort_values("ADAR", ascending=False)
        sig_matrix = sig_matrix.loc[matrix.index]

    annotations = pd.DataFrame(index=matrix.index, columns=matrix.columns, dtype=str)
    for cancer_type in matrix.index:
        for gene in matrix.columns:
            r_value = matrix.loc[cancer_type, gene]
            stars = sig_matrix.loc[cancer_type, gene]
            stars = "" if stars == "ns" else stars
            annotations.loc[cancer_type, gene] = f"{r_value:.2f}{stars}"

    height = max(6, 0.28 * len(matrix))
    fig, ax = plt.subplots(figsize=(6, height))
    sns.heatmap(
        matrix,
        annot=annotations,
        fmt="",
        cmap="RdBu_r",
        center=0,
        vmin=-0.5,
        vmax=0.8,
        linewidths=0.5,
        cbar_kws={"label": "Pearson r"},
        ax=ax,
    )
    ax.set_title("AEI vs ADAR family expression (per cancer type)", fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("Cancer type")
    ax.text(
        1.52, 1.02,
        "Pearson correlation\n* p < 0.05\n** p < 0.01\n*** p < 0.001",
        transform=ax.transAxes,
        ha="left", va="top", fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85, edgecolor="0.8"),
    )
    plt.tight_layout()
    save_figure(SCRIPT, os.path.join(outdir, "AEI_vs_ADAR_family_heatmap.png"), bbox_inches="tight")


def figure3_tumor_vs_normal(tumor, normal, outdir):
    """Side-by-side boxplots of AEI in tumor and batch-corrected matched normal samples."""
    combined, batch_info = get_matched_tumor_normal_data(
        tumor, normal, batch_correct=True,
    )
    if combined.empty:
        print(f"[{SCRIPT}] Skipped: tumor vs normal figure (no matched normal reference)")
        return

    if batch_info is not None:
        print(f"[{SCRIPT}] Applied AEI batch correction: global median="
              f"{batch_info['global_median']:.3f}, tumor {batch_info['tumor_shift']:+.3f}, "
              f"normal {batch_info['normal_shift']:+.3f}")

    combined = clean_labels(combined)
    order = sorted(combined["Label"].unique())

    fig, ax = plt.subplots(figsize=(12, 5))
    sns.boxplot(
        data=combined, x="Label", y="AEI", hue="Status", order=order,
        palette={"tumor": "salmon", "normal": "lightblue"}, ax=ax,
        flierprops=dict(marker="o", markersize=2, alpha=0.3),
    )

    y_range = combined["AEI"].max() - combined["AEI"].min()
    y_offset = 0.03 * y_range if y_range > 0 else 0.05
    y_top = combined.groupby("Label")["AEI"].max()

    for index, label in enumerate(order):
        subset = combined[combined["Label"] == label]
        tumor_aei = subset[subset["Status"] == "tumor"]["AEI"]
        normal_aei = subset[subset["Status"] == "normal"]["AEI"]
        if len(tumor_aei) == 0 or len(normal_aei) == 0:
            continue
        _, p_value = stats.mannwhitneyu(tumor_aei, normal_aei, alternative="two-sided")
        ax.text(
            index,
            y_top[label] + y_offset,
            significance_stars(p_value),
            ha="center",
            va="bottom",
            fontsize=11,
        )

    plt.xticks(rotation=90)
    ax.set_title("AEI in tumor vs matched normal tissue (batch-corrected)", fontsize=13)
    ax.set_xlabel("Cancer type", fontsize=11)
    ax.set_ylabel("AEI (%)", fontsize=11)
    ax.legend(title="Status", loc="upper left")
    add_significance_legend(ax, loc="upper right")
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

    correlation_tables = _load_correlation_tables(args.correlations)
    for gene in EXPRESSION_GENES:
        if gene not in tumor.columns:
            print(f"[{SCRIPT}] Skipped: AEI vs {gene} figure ({gene} column missing)")
            continue
        if gene not in correlation_tables:
            print(f"[{SCRIPT}] Skipped: AEI vs {gene} figure (AEI_vs_{gene}.csv missing)")
            continue
        figure_aei_vs_expression(tumor, correlation_tables[gene], gene, args.output)

    figure_aei_expression_heatmap(correlation_tables, args.output)

    figure3_tumor_vs_normal(tumor, normal, args.output)

    print(f"[{SCRIPT}] Done. Output directory: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
