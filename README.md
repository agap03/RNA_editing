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
│   ├── AEI_analysis.py
│   ├── REDI_statistics_figures.py
│   ├── combined_analysis_figures.py
│   ├── survival_analysis.py
│   ├── survival_analysis_summary.py
│   ├── redi_utils.py         # Shared REDI data loading
│   └── plot_utils.py         # Shared figure export
├── data/                     # Input files
└── outputs/                  # Generated results
```

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

## Data setup

Input data is present in the data folder.

## Running the pipeline

From the project root:

```bash
conda activate projekt
snakemake --cores 
```

## Workflow stages

1. **Context extraction** — Annotate each editing site with trinucleotide context (COSMIC notation) from `hg38.fa`.
2. **Trinucleotide analysis** — Filter low-confidence sites (bottom 10% by `edSample`), build 192-context spectra.
3. **Signature fitting** — Decompose spectra onto COSMIC RNA-SBS signatures via SigProfilerAssignment.
4. **AEI analysis** — Cross-cancer comparisons, tumour vs normal tests, AEI vs *ADAR*/*ADARB1*/*ADARB2* correlations.
5. **REDI figures** — Boxplots, scatter plots, tumour–normal comparisons.
6. **Combined figures** — RNA-SBS1 activity vs mean AEI and vs AEI–ADAR correlation.
7. **Survival analysis** — Per-cancer Kaplan–Meier curves (high vs low AEI, median split) and log-rank tests.

## Outputs

| Directory | Contents |
|-----------|----------|
| `outputs/context_extraction_output/` | `{CANCER}_contexts.csv` |
| `outputs/trinucleotide_analysis_output/` | `{CANCER}_spectrum.csv`, profile PNGs |
| `outputs/signature_fitting_output/` | Signature activities, decomposition stats |
| `outputs/AEI_analysis_output/` | Statistical CSV tables |
| `outputs/REDI_statistics_figures_output/` | AEI and expression figures |
| `outputs/combined_analysis_figures_output/` | Integration scatter plots |
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