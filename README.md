# RNA Editing Analysis Pipeline

A Snakemake workflow integrating **REIA** RNA editing mutation spectra, **REDIportal** Adenosine Editing Index (AEI) statistics, and **TCGA** overall survival data across 33 cancer types.

## Overview

The pipeline has a few separate parts:

| Part | Input | Analysis |
|--------|-------|----------|
| **REIA** | Per-cancer editing mutation CSVs | Trinucleotide context → spectrum → COSMIC RNA-SBS signature fitting |
| **REDI** | REDIportal tumour/normal tables | AEI statistics, ADAR family correlations, figures |
| **Clinical** | TCGA survival TSVs | Kaplan–Meier survival by high vs low AEI |
| **Integration** | Signature activities + AEI | RNA-SBS1 activity vs mean AEI and AEI–ADAR correlation |

## Project structure

```
projekt/
├── Snakefile                 # Workflow definition
├── requirements.txt          # Python dependencies
├── scripts/
│   ├── context_extraction.py
│   ├── trinucleotide_analysis.py
│   ├── signature_fitting.py
│   ├── signature_figures.py
│   ├── AEI_analysis.py
│   ├── REDI_statistics_figures.py
│   ├── combined_analysis_figures.py
│   ├── survival_analysis.py
│   ├── survival_analysis_summary.py
│   ├── redi_utils.py         # Shared REDI data loading
│   └── plot_utils.py         # Shared figure export
├── data/
│   ├── REIA_data/
│   │   └── {CANCER}.csv                    # 33 per-cancer editing mutation tables
│   ├── REDI_data/
│   │   ├── cancer_REDIportal.csv
│   │   └── normal_REDIportal.csv
│   ├── clinical_data/
│   │   └── TCGA-{CANCER}.survival.tsv      # 32 TCGA survival files
│   ├── hg38.fa                               # Reference genome (download separately)
│   └── COSMIC_Human_RNA-SBS-192_GRCh37_v3.6.csv
└── outputs/
    ├── context_extraction_output/
    │   └── {CANCER}_contexts.csv
    ├── trinucleotide_analysis_output/
    │   ├── {CANCER}_spectrum.csv
    │   └── {CANCER}_trinucleotide_profile.png
    ├── signature_fitting_output/
    │   ├── JOB_METADATA_SPA.txt
    │   └── Decompose_Solution/
    │       ├── Activities/
    │       │   ├── Decompose_Solution_Activities.txt
    │       │   └── Decomposed_MutationType_Probabilities.txt
    │       ├── Signatures/
    │       │   └── Decompose_Solution_Signatures.txt
    │       ├── Solution_Stats/
    │       │   ├── Decompose_Solution_Samples_Stats.txt
    │       │   ├── Decompose_Solution_Signature_Assignment_log.txt
    │       │   └── Cosmic_SBS192_Decomposition_Log.txt
    │       └── De_Novo_map_to_COSMIC_SBS192.csv
    ├── signature_figures_output/
    │   └── activity_by_cancer_type.png
    ├── AEI_analysis_output/
    │   ├── AEI_summary_statistics.csv
    │   ├── AEI_pairwise.csv
    │   ├── AEI_tumor_vs_normal.csv
    │   ├── AEI_vs_ADAR.csv
    │   ├── AEI_vs_ADARB1.csv
    │   └── AEI_vs_ADARB2.csv
    ├── REDI_statistics_figures_output/
    │   ├── AEI_across_cancers.png
    │   ├── AEI_vs_ADAR.png
    │   ├── AEI_vs_ADARB1.png
    │   ├── AEI_vs_ADARB2.png
    │   ├── AEI_vs_ADAR_family_heatmap.png
    │   └── tumor_vs_normal.png
    ├── combined_analysis_figures_output/
    │   ├── activity_vs_mean_AEI.png
    │   └── activity_vs_pearson_correlation.png
    └── survival_analysis_output/
        ├── {CANCER}_merged.csv
        ├── {CANCER}_km.png
        ├── {CANCER}_logrank.csv
        └── survival_logrank_results.csv
```

`{CANCER}` denotes one cancer type per file (e.g. `BLCA`, `LUAD`, `STAD`).

## Requirements

- Python 3.11+
- Snakemake
- Linux or WSL2

Install dependencies:

```bash
conda create -n env python=3.11 -y
conda activate env
pip install -r requirements.txt
```

Download and decompress the human reference genome (required for context extraction):

```bash
mkdir -p data
curl -L -o data/hg38.fa.gz https://hgdownload.gi.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip data/hg38.fa.gz
```

This produces `data/hg38.fa` (~3 GB uncompressed).

## Data setup

Other input files should be placed in `data/` (REIA, REDIportal, TCGA clinical, COSMIC signatures). See `.gitignore` for paths not tracked in git.

## Running the pipeline

From the project root:

```bash
conda activate env
snakemake --cores 
```

## Workflow stages

1. **Context extraction** — Annotate each editing site with trinucleotide context (COSMIC notation) from `hg38.fa`.
2. **Trinucleotide analysis** — Filter low-confidence sites (bottom 10% by `edSample`), build 192-context spectra.
3. **Signature fitting** — Decompose spectra onto COSMIC RNA-SBS signatures via SigProfilerAssignment.
4. **Signature figures** — Summarise RNA-SBS1 activity across cancer types.
5. **AEI analysis** — Cross-cancer comparisons, tumour vs normal tests, AEI vs *ADAR*/*ADARB1*/*ADARB2* correlations.
6. **REDI figures** — Boxplots, per-gene scatter plots, ADAR-family correlation heatmap, tumour–normal comparisons.
7. **Combined figures** — Integrate REIA signature activity with REDI AEI metrics (two scatter plots).
8. **Survival analysis** — Per-cancer Kaplan–Meier curves (high vs low AEI, median split) and log-rank tests.

## Outputs

| Directory | Contents |
|-----------|----------|
| `outputs/context_extraction_output/` | `{CANCER}_contexts.csv` |
| `outputs/trinucleotide_analysis_output/` | `{CANCER}_spectrum.csv`, profile PNGs |
| `outputs/signature_fitting_output/` | Signature activities, decomposition stats |
| `outputs/signature_figures_output/` | `activity_by_cancer_type.png` |
| `outputs/AEI_analysis_output/` | Statistical CSV tables |
| `outputs/REDI_statistics_figures_output/` | AEI boxplots, expression scatter plots, `AEI_vs_ADAR_family_heatmap.png`, tumour–normal figure |
| `outputs/combined_analysis_figures_output/` | `activity_vs_mean_AEI.png`, `activity_vs_pearson_correlation.png` |
| `outputs/survival_analysis_output/` | Per-cancer KM plots, merged data, log-rank summary |

## Running individual scripts

Each script can be run standalone:

```bash
python scripts/AEI_analysis.py \
    --cancers data/REDI_data/cancer_REDIportal.csv \
    --normal data/REDI_data/normal_REDIportal.csv \
    --output outputs/AEI_analysis_output
```

```bash
python scripts/survival_analysis.py \
    --cancers data/REDI_data/cancer_REDIportal.csv \
    --clinical data/clinical_data/TCGA-STAD.survival.tsv \
    --output outputs/survival_analysis_output
```
