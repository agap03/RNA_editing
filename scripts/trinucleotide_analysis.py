#!/usr/bin/env python3
"""
Create a trinucleotide mutation spectrum from one contexts CSV.

Filters low-confidence editing sites, counts Context occurrences, and saves
both the spectrum table and a diagnostic bar plot for the cancer sample.
"""
import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

from plot_utils import save_figure

SCRIPT = "trinucleotide_analysis"


def plot_spectrum(df, path):
    """Horizontal bar plot of context counts, sorted for readability."""
    plt.figure(figsize=(6, 5))
    df_sorted = df.sort_values("Count")
    plt.barh(df_sorted["Context"], df_sorted["Count"])
    plt.xlabel("Count")
    plt.ylabel("Context")
    plt.tight_layout()
    save_figure(SCRIPT, path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input CSV file with Context column")
    parser.add_argument("-o", "--output", required=True, help="Output spectrum CSV path")
    parser.add_argument("-p", "--plot", required=True, help="Output trinucleotide profile PNG path")
    args = parser.parse_args()

    sample = os.path.splitext(os.path.basename(args.input))[0].replace("_contexts", "")
    print(f"[{SCRIPT}] Sample: {sample}")

    mut = pd.read_csv(args.input)

    # Remove the bottom 10% of sites by edSample to reduce noise from weak editing calls.
    number_threshold = mut["edSample"].quantile(0.10)
    mut = mut[mut["edSample"] > number_threshold]

    if "Context" not in mut.columns:
        raise ValueError(f"No Context column in {args.input}")

    spectrum = mut["Context"].value_counts().reset_index()
    spectrum.columns = ["Context", "Count"]

    for path in (args.output, args.plot):
        out_dir = os.path.dirname(path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    spectrum.to_csv(args.output, index=False)
    plot_spectrum(spectrum, args.plot)
    print(f"[{SCRIPT}] Saved: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
