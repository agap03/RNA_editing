#!/usr/bin/env python3
"""
Fit COSMIC RNA-SBS-192 signatures to per-cancer trinucleotide spectra.

Collects all *_spectrum.csv files, aligns them to the full 192-context space,
and runs SigProfilerAssignment decompose_fit to estimate signature activities.
"""

import argparse
import glob
import os

import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze

SCRIPT = "signature_fitting"


def load_signatures(sig_file):
    """Load the reference signature matrix with contexts as row labels."""
    sig = pd.read_csv(sig_file)
    if "Type" not in sig.columns:
        raise ValueError("Signature matrix must contain a 'Type' column")
    sig = sig.set_index("Type")
    return sig.astype(float)


def load_spectra(input_folder, contexts):
    """
    Load per-cancer spectra and force each sample into the same 192-context
    space as the reference signatures.

    Missing contexts are filled with zero so SigProfiler receives complete vectors.
    """
    all_files = glob.glob(os.path.join(input_folder, "*_spectrum.csv"))
    if not all_files:
        raise ValueError(f"No *_spectrum.csv files found in {input_folder}")

    sample_matrices = []
    for f in all_files:
        sample_name = os.path.basename(f).replace("_spectrum.csv", "")
        df = pd.read_csv(f)

        if "Context" not in df.columns or "Count" not in df.columns:
            raise ValueError(f"{f} must contain Context and Count columns")

        df = df.set_index("Context").reindex(contexts).fillna(0)
        df = df.rename(columns={"Count": sample_name})
        sample_matrices.append(df)

    return pd.concat(sample_matrices, axis=1).fillna(0)


def run_assignment(samples, signatures, output_dir):
    """Run signature decomposition and write SigProfiler outputs."""
    os.makedirs(output_dir, exist_ok=True)

    signatures = signatures.reindex(index=samples.index).fillna(0)
    if samples.shape[0] != signatures.shape[0]:
        raise ValueError(
            f"Dimension mismatch!\n"
            f"samples: {samples.shape}\n"
            f"signatures: {signatures.shape}"
        )

    return Analyze.decompose_fit(
        samples=samples,
        signatures=signatures,
        output=output_dir,
        verbose=True,
        collapse_to_SBS96=False,  # Keep full RNA-SBS-192 resolution.
    )


def main():
    parser = argparse.ArgumentParser(description="SigProfilerAssignment RNA-SBS-192 fitting")
    parser.add_argument("-i", "--input", required=True, help="Input folder with *_spectrum.csv files")
    parser.add_argument("-o", "--output", required=True, help="Output folder")
    parser.add_argument("-s", "--signatures", required=True, help="RNA-SBS-192 signature matrix CSV")
    args = parser.parse_args()

    print(f"[{SCRIPT}] Loading signatures...")
    signatures = load_signatures(args.signatures)
    print(f"[{SCRIPT}] Signature contexts: {len(signatures.index)}")

    print(f"[{SCRIPT}] Loading spectra...")
    samples = load_spectra(args.input, signatures.index.tolist())
    print(f"[{SCRIPT}] Loaded {samples.shape[1]} samples ({samples.shape[0]} contexts)")

    print(f"[{SCRIPT}] Running SigProfilerAssignment...")
    run_assignment(samples, signatures, args.output)

    print(f"[{SCRIPT}] Done. Output directory: {os.path.abspath(args.output)}")


if __name__ == "__main__":
    main()
