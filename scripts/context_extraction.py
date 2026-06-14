#!/usr/bin/env python3
"""
Extract trinucleotide sequence context for point mutations.

Reads a per-cancer REIA CSV with columns Chr, Position(1base), Ref, Ed, Strand
and writes the same rows with added Trinucleotide and Context columns in
COSMIC RNA-SBS-192 notation, e.g. A[C>T]G.
"""
import argparse
import os

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

SCRIPT = "context_extraction"


def rev_comp_base(b):
    return str(Seq(b).reverse_complement())


def rev_comp_context(s):
    return str(Seq(s).reverse_complement())


def get_contexts(df, seq_index):
    """
    Look up the 3 bp reference sequence around each mutation and build
    a COSMIC-style context string.

    Negative-strand variants are reverse-complemented so all contexts are
    reported on the transcriptional forward strand.
    """
    tri_list = []
    contexts = []

    for _, row in df.iterrows():
        chrom = str(row["Chr"])
        pos = int(row["Position(1base)"])

        # Mutation position is 1-based; extract [pos-2, pos+1) from the genome.
        start = pos - 2
        end = pos + 1

        tri = ""

        # FASTA records may use either "chr1" or "1" naming.
        candidates = [chrom]
        if chrom.startswith("chr"):
            candidates.append(chrom.replace("chr", "", 1))
        else:
            candidates.append("chr" + chrom)

        for c in candidates:
            try:
                seq = seq_index[c].seq
            except KeyError:
                continue

            s = str(seq[max(0, start):end]).upper()
            if len(s) == 3:
                tri = s
                break

        if len(tri) != 3:
            tri_list.append(None)
            contexts.append(None)
            continue

        ref = str(row["Ref"]).upper()
        alt = str(row["Ed"]).upper()
        strand = str(row.get("Strand", "+"))

        context = tri
        if strand == "-":
            context = rev_comp_context(context)
            ref = rev_comp_base(ref)
            alt = rev_comp_base(alt)

        ctxt = f"{context[0]}[{ref}>{alt}]{context[2]}"
        tri_list.append(context)
        contexts.append(ctxt)

    df["Trinucleotide"] = tri_list
    df["Context"] = contexts
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input CSV file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    parser.add_argument("-g", "--genome", required=True, help="Indexed FASTA (matching chromosome names)")
    args = parser.parse_args()

    # Indexed FASTA avoids loading the whole genome into memory.
    seq_index = SeqIO.index(args.genome, "fasta")

    sample = os.path.splitext(os.path.basename(args.input))[0]
    print(f"[{SCRIPT}] Sample: {sample}")

    df = pd.read_csv(args.input)
    if not {"Chr", "Position(1base)", "Ref", "Ed"}.issubset(df.columns):
        raise ValueError(f"Missing required columns in {args.input}")

    out = get_contexts(df, seq_index)

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f"[{SCRIPT}] Saved: {os.path.abspath(args.output)}")
    out.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
