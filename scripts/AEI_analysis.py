#!/usr/bin/env python3
"""
Run statistical analyses on REDIportal tumor and normal samples.

Produces summary tables and per-cancer correlation statistics used downstream
by the REDI figure scripts and combined analysis figures.
"""

import argparse
import os

import pandas as pd
from itertools import combinations
from scipy import stats
from statsmodels.stats.multitest import multipletests

from plot_utils import significance_stars
from redi_utils import NORMAL_TISSUE_MAP, get_matched_tumor_normal_data, load_redi_data

SCRIPT = "AEI_analysis"

# Minimum samples required for per-cancer tests.
MIN_SAMPLES = 5


def analysis1_aei_across_cancers(tumor):
    """Summarise AEI by cancer type and run pairwise Mann-Whitney U tests."""
    print(f"[{SCRIPT}] Analysis 1: AEI across cancer types")

    summary = tumor.groupby("Cancer_type")["AEI"].agg(["median", "mean", "std", "count"]).round(3)
    print(f"[{SCRIPT}] Summary statistics:")
    print(summary.to_string())

    cancer_types = tumor["Cancer_type"].unique()
    pairs = list(combinations(cancer_types, 2))

    valid_pairs = []
    p_values = []
    pair_meta = []
    for ct1, ct2 in pairs:
        g1 = tumor[tumor["Cancer_type"] == ct1]["AEI"].dropna()
        g2 = tumor[tumor["Cancer_type"] == ct2]["AEI"].dropna()
        if len(g1) < MIN_SAMPLES or len(g2) < MIN_SAMPLES:
            pair_meta.append((ct1, ct2, None))
            continue
        _, p_pair = stats.mannwhitneyu(g1, g2, alternative="two-sided")
        valid_pairs.append((ct1, ct2))
        p_values.append(p_pair)
        pair_meta.append((ct1, ct2, p_pair))

    # Benjamini-Hochberg correction across all valid pairwise comparisons.
    if len(p_values) > 0:
        _, p_corr, _, _ = multipletests(p_values, method="fdr_bh")
        p_corr_iter = iter(p_corr)
    else:
        p_corr_iter = iter([])

    print(f"[{SCRIPT}] Pairwise comparisons (BH-corrected):")
    pairwise_results = []
    for meta in pair_meta:
        ct1, ct2, raw_p = meta
        if raw_p is None:
            print(f"[{SCRIPT}]   {ct1} vs {ct2}: skipped (insufficient samples)")
            pairwise_results.append({
                "Comparison": f"{ct1}_vs_{ct2}",
                "p_adjusted": None,
                "Significance": "insufficient_samples",
            })
            continue
        p_adj = next(p_corr_iter)
        sig = significance_stars(p_adj)
        print(f"[{SCRIPT}]   {ct1} vs {ct2}: p={p_adj:.4f} {sig}")
        pairwise_results.append({
            "Comparison": f"{ct1}_vs_{ct2}",
            "p_adjusted": p_adj,
            "Significance": sig,
        })

    return summary, pd.DataFrame(pairwise_results)


def analysis_aei_vs_expression(tumor, gene):
    """Pearson correlation between AEI and an editing-gene expression column."""
    print(f"[{SCRIPT}] AEI vs {gene} expression correlation")

    correlation_results = []
    for ct in sorted(tumor["Cancer_type"].unique()):
        subset = tumor[tumor["Cancer_type"] == ct].dropna(subset=["AEI", gene])
        if len(subset) < MIN_SAMPLES:
            print(f"[{SCRIPT}]   {ct}: skipped (insufficient samples)")
            continue
        r, p = stats.pearsonr(subset[gene], subset["AEI"])
        sig = significance_stars(p)
        print(f"[{SCRIPT}]   {ct}: r={r:.3f}, p={p:.4f} {sig}, n={len(subset)}")
        correlation_results.append({
            "Cancer_type": ct,
            "Pearson_r": round(r, 3),
            "p_value": round(p, 6),
            "Significance": sig,
            "n": len(subset),
        })

    return pd.DataFrame(correlation_results)


def analysis3_tumor_vs_normal(tumor, normal, normal_map):
    """Compare AEI between tumor and batch-corrected matched normal tissue."""
    print(f"[{SCRIPT}] Analysis 3: Tumor vs normal AEI")

    combined, batch_info = get_matched_tumor_normal_data(
        tumor, normal, normal_map, batch_correct=True,
    )
    if combined.empty:
        print(f"[{SCRIPT}]   skipped (no matched normal reference)")
        return pd.DataFrame()

    if batch_info is not None:
        print(f"[{SCRIPT}]   batch correction: global median={batch_info['global_median']:.3f}, "
              f"tumor shift={batch_info['tumor_shift']:+.3f}, "
              f"normal shift={batch_info['normal_shift']:+.3f}")

    tumor_normal_results = []
    for cancer_type in sorted(combined["Cancer_type"].unique()):
        subset = combined[combined["Cancer_type"] == cancer_type]
        tumor_aei = subset[subset["Status"] == "tumor"]["AEI"]
        normal_aei = subset[subset["Status"] == "normal"]["AEI"]

        _, p = stats.mannwhitneyu(tumor_aei, normal_aei, alternative="two-sided")
        direction = "higher" if tumor_aei.median() > normal_aei.median() else "lower"
        sig = significance_stars(p)

        print(f"[{SCRIPT}]   {cancer_type}: tumor={tumor_aei.median():.3f}, "
              f"normal={normal_aei.median():.3f}, "
              f"{direction} in tumor, p={p:.4f} {sig}")

        tumor_normal_results.append({
            "Cancer_type": cancer_type,
            "Tumor_median_AEI": round(tumor_aei.median(), 3),
            "Normal_median_AEI": round(normal_aei.median(), 3),
            "Direction": direction,
            "p_value": round(p, 6),
            "Significance": sig,
            "n_tumor": len(tumor_aei),
            "n_normal": len(normal_aei),
        })

    return pd.DataFrame(tumor_normal_results)


def main():
    parser = argparse.ArgumentParser(description="Run statistical analyses on REDIportal data")
    parser.add_argument("--cancers", "-c", required=True, help="Path to cancer_REDIportal.csv")
    parser.add_argument("--normal", "-n", required=True, help="Path to normal_REDIportal.csv")
    parser.add_argument("--output", "-o", required=True, help="Output directory for analysis results")
    args = parser.parse_args()

    print(f"[{SCRIPT}] Loading cancer and normal data...")
    tumor, normal = load_redi_data(args.cancers, args.normal, script=SCRIPT)

    os.makedirs(args.output, exist_ok=True)

    summary_df, pairwise_df = analysis1_aei_across_cancers(tumor)
    tumor_normal_df = analysis3_tumor_vs_normal(tumor, normal, NORMAL_TISSUE_MAP)

    output_files = {
        "AEI_summary_statistics.csv": summary_df,
        "AEI_pairwise.csv": pairwise_df,
        "AEI_tumor_vs_normal.csv": tumor_normal_df,
    }
    for gene in ["ADAR", "ADARB1", "ADARB2"]:
        if gene not in tumor.columns:
            print(f"[{SCRIPT}] Skipped: AEI vs {gene} ({gene} column missing)")
            continue
        output_files[f"AEI_vs_{gene}.csv"] = analysis_aei_vs_expression(tumor, gene)

    for filename, df in output_files.items():
        path = os.path.join(args.output, filename)
        if filename == "AEI_summary_statistics.csv":
            df.to_csv(path)
        else:
            df.to_csv(path, index=False)
        print(f"[{SCRIPT}] Saved: {os.path.abspath(path)}")

    print(f"[{SCRIPT}] Done. Output directory: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
