#!/usr/bin/env python3
"""
Combine per-cancer survival log-rank results into one summary table.

This script is called once by Snakemake after all individual survival_analysis
jobs have finished.
"""
import argparse
import glob
import os

import pandas as pd

SCRIPT = "survival_analysis_summary"


def main():
    parser = argparse.ArgumentParser(description="Summarize per-cancer survival results")
    parser.add_argument("--input", "-i", required=True, help="Directory with *_logrank.csv files")
    parser.add_argument("--output", "-o", required=True, help="Output summary CSV path")
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.input, "*_logrank.csv")))
    if not files:
        raise ValueError(f"No *_logrank.csv files found in {args.input}")

    summary = pd.concat((pd.read_csv(path) for path in files), ignore_index=True)
    summary = summary.sort_values("Cancer_type").reset_index(drop=True)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    summary.to_csv(args.output, index=False)
    print(f"[{SCRIPT}] Saved: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
