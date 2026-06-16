#!/usr/bin/env python3
"""
Generate summary figures from COSMIC RNA-SBS signature fitting.

This script belongs to the REIA branch and uses only signature-fitting output.
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

from plot_utils import save_figure

SCRIPT = "signature_figures"


def load_rna_sbs1_activity(activity_path):
    """Load per-cancer RNA-SBS1 activity from SigProfiler activities file."""
    activity = pd.read_csv(activity_path, sep="\t")
    activity["Cancer_type"] = activity["Samples"]
    return activity[["Cancer_type", "0"]].rename(columns={"0": "Activity"})


def figure_activity_by_cancer_type(activity, outdir):
    """Bar chart of RNA-SBS1 signature activity across cancer types."""
    activity = activity.sort_values("Cancer_type")

    plt.figure(figsize=(8, 6))
    plt.bar(activity["Cancer_type"], activity["Activity"])
    plt.xticks(rotation=90)
    plt.xlabel("Cancer type")
    plt.ylabel("Activity")
    plt.title("Activity of signature RNA-SBS1 across cancer types")
    plt.tight_layout()
    save_figure(SCRIPT, os.path.join(outdir, "activity_by_cancer_type.png"))


def main():
    parser = argparse.ArgumentParser(description="Generate figures from signature fitting output")
    parser.add_argument("--activity", "-a", required=True, help="Path to Decompose_Solution_Activities.txt")
    parser.add_argument("--output", "-o", required=True, help="Output directory for figures")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    print(f"[{SCRIPT}] Loading signature activities...")
    activity = load_rna_sbs1_activity(args.activity)
    figure_activity_by_cancer_type(activity, args.output)

    print(f"[{SCRIPT}] Done. Output directory: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
