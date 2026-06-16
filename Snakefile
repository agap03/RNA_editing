# Combined workflow:
# 1. Context extraction
# 2. Trinucleotide spectrum generation
# 3. Signature fitting
# 4. Signature figures 
# 5. AEI analyses
# 6. REDI statistics figures
# 7. Combined analysis figures 
# 8. Survival analysis

SAMPLES = glob_wildcards("data/REIA_data/{sample}.csv").sample
CLINICAL = glob_wildcards("data/clinical_data/TCGA-{cancer}.survival.tsv").cancer

GENOME = "data/hg38.fa"
SIGNATURES = "data/COSMIC_Human_RNA-SBS-192_GRCh37_v3.6.csv"


rule all:
    input:
        # Context extraction
        expand(
            "outputs/context_extraction_output/{sample}_contexts.csv",
            sample=SAMPLES
        ),

        # Trinucleotide analysis
        expand(
            "outputs/trinucleotide_analysis_output/{sample}_spectrum.csv",
            sample=SAMPLES
        ),

        # Signature fitting
        "outputs/signature_fitting_output/Decompose_Solution/Activities/Decompose_Solution_Activities.txt",

        # Signature figures
        "outputs/signature_figures_output/activity_by_cancer_type.png",

        # AEI analysis outputs
        "outputs/AEI_analysis_output/AEI_pairwise.csv",
        "outputs/AEI_analysis_output/AEI_summary_statistics.csv",
        "outputs/AEI_analysis_output/AEI_vs_ADAR.csv",
        "outputs/AEI_analysis_output/AEI_vs_ADARB1.csv",
        "outputs/AEI_analysis_output/AEI_vs_ADARB2.csv",
        "outputs/AEI_analysis_output/AEI_tumor_vs_normal.csv",

        # REDI statistics figures
        "outputs/REDI_statistics_figures_output/AEI_across_cancers.png",
        "outputs/REDI_statistics_figures_output/AEI_vs_ADAR.png",
        "outputs/REDI_statistics_figures_output/AEI_vs_ADARB1.png",
        "outputs/REDI_statistics_figures_output/AEI_vs_ADARB2.png",
        "outputs/REDI_statistics_figures_output/AEI_vs_ADAR_family_heatmap.png",
        "outputs/REDI_statistics_figures_output/tumor_vs_normal.png",

        # Combined analysis figures
        "outputs/combined_analysis_figures_output/activity_vs_mean_AEI.png",
        "outputs/combined_analysis_figures_output/activity_vs_pearson_correlation.png",

        # Survival analysis
        expand(
            "outputs/survival_analysis_output/{cancer}_km.png",
            cancer=CLINICAL
        ),
        "outputs/survival_analysis_output/survival_logrank_results.csv",


# Context extraction

rule context_extraction:
    input:
        "data/REIA_data/{sample}.csv"
    output:
        "outputs/context_extraction_output/{sample}_contexts.csv"
    params:
        genome=GENOME
    shell:
        """
        python scripts/context_extraction.py \
            --input {input} \
            --output {output} \
            --genome {params.genome}
        """


# Trinucleotide analysis

rule trinucleotide_analysis:
    input:
        "outputs/context_extraction_output/{sample}_contexts.csv"
    output:
        spectrum="outputs/trinucleotide_analysis_output/{sample}_spectrum.csv",
        plot="outputs/trinucleotide_analysis_output/{sample}_trinucleotide_profile.png"
    shell:
        """
        python scripts/trinucleotide_analysis.py \
            --input {input} \
            --output {output.spectrum} \
            --plot {output.plot}
        """


# Signature fitting

rule signature_fitting:
    input:
        spectra=expand(
            "outputs/trinucleotide_analysis_output/{sample}_spectrum.csv",
            sample=SAMPLES
        ),
        signatures=SIGNATURES
    output:
        activities="outputs/signature_fitting_output/Decompose_Solution/Activities/Decompose_Solution_Activities.txt",
        metadata="outputs/signature_fitting_output/JOB_METADATA_SPA.txt"
    shell:
        """
        python scripts/signature_fitting.py \
            --input outputs/trinucleotide_analysis_output \
            --output outputs/signature_fitting_output \
            --signatures {input.signatures}
        """


# Signature figures 

rule signature_figures:
    input:
        activity="outputs/signature_fitting_output/Decompose_Solution/Activities/Decompose_Solution_Activities.txt",
    output:
        activity_by_cancer="outputs/signature_figures_output/activity_by_cancer_type.png",
    shell:
        """
        python scripts/signature_figures.py \
            --activity {input.activity} \
            --output outputs/signature_figures_output
        """


# AEI analyses

rule aei_analysis:
    input:
        cancers="data/REDI_data/cancer_REDIportal.csv",
        normal="data/REDI_data/normal_REDIportal.csv"
    output:
        pairwise="outputs/AEI_analysis_output/AEI_pairwise.csv",
        summary="outputs/AEI_analysis_output/AEI_summary_statistics.csv",
        correlation_adar="outputs/AEI_analysis_output/AEI_vs_ADAR.csv",
        correlation_adarb1="outputs/AEI_analysis_output/AEI_vs_ADARB1.csv",
        correlation_adarb2="outputs/AEI_analysis_output/AEI_vs_ADARB2.csv",
        tumor_vs_normal="outputs/AEI_analysis_output/AEI_tumor_vs_normal.csv"
    shell:
        """
        python scripts/AEI_analysis.py \
            --cancers {input.cancers} \
            --normal {input.normal} \
            --output outputs/AEI_analysis_output
        """


# REDI statistics figures

rule redi_statistics_figures:
    input:
        cancers="data/REDI_data/cancer_REDIportal.csv",
        normal="data/REDI_data/normal_REDIportal.csv",
        correlations=rules.aei_analysis.output,
    output:
        across_cancers="outputs/REDI_statistics_figures_output/AEI_across_cancers.png",
        aei_vs_adar="outputs/REDI_statistics_figures_output/AEI_vs_ADAR.png",
        aei_vs_adarb1="outputs/REDI_statistics_figures_output/AEI_vs_ADARB1.png",
        aei_vs_adarb2="outputs/REDI_statistics_figures_output/AEI_vs_ADARB2.png",
        aei_heatmap="outputs/REDI_statistics_figures_output/AEI_vs_ADAR_family_heatmap.png",
        tumor_vs_normal="outputs/REDI_statistics_figures_output/tumor_vs_normal.png"
    shell:
        """
        python scripts/REDI_statistics_figures.py \
            --cancers {input.cancers} \
            --normal {input.normal} \
            --correlations outputs/AEI_analysis_output \
            --output outputs/REDI_statistics_figures_output
        """


# Combined analysis figures

rule combined_analysis_figures:
    input:
        activity="outputs/signature_fitting_output/Decompose_Solution/Activities/Decompose_Solution_Activities.txt",
        correlation="outputs/AEI_analysis_output/AEI_vs_ADAR.csv",
        aei="outputs/AEI_analysis_output/AEI_summary_statistics.csv"
    output:
        activity_vs_aei="outputs/combined_analysis_figures_output/activity_vs_mean_AEI.png",
        activity_vs_correlation="outputs/combined_analysis_figures_output/activity_vs_pearson_correlation.png"
    shell:
        """
        python scripts/combined_analysis_figures.py \
            --activity {input.activity} \
            --correlation {input.correlation} \
            --aei {input.aei} \
            --output outputs/combined_analysis_figures_output
        """


# Survival analysis

rule survival_analysis:
    input:
        cancers="data/REDI_data/cancer_REDIportal.csv",
        clinical="data/clinical_data/TCGA-{cancer}.survival.tsv",
    output:
        merged="outputs/survival_analysis_output/{cancer}_merged.csv",
        km_plot="outputs/survival_analysis_output/{cancer}_km.png",
        logrank="outputs/survival_analysis_output/{cancer}_logrank.csv",
    shell:
        """
        python scripts/survival_analysis.py \
            --cancers {input.cancers} \
            --clinical {input.clinical} \
            --output outputs/survival_analysis_output
        """


rule survival_analysis_summary:
    input:
        expand(
            "outputs/survival_analysis_output/{cancer}_logrank.csv",
            cancer=CLINICAL
        ),
    output:
        "outputs/survival_analysis_output/survival_logrank_results.csv",
    shell:
        """
        python scripts/survival_analysis_summary.py \
            --input outputs/survival_analysis_output \
            --output {output}
        """
