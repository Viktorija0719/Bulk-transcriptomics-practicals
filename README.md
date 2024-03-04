# Bulk-transcriptomics-practicals
## RNA-seq Analysis Workflow
###Project Description
This project houses a Snakemake workflow designed for comprehensive RNA sequencing (RNA-seq) data analysis. The workflow encompasses several key steps: quality control with FastQC and MultiQC, adapter trimming using bbduk, and read alignment with the STAR aligner. 
Workflow Steps:
    1. Quality Evaluation: Assess the quality of raw sequencing reads using FastQC.
    2. Quality Report Aggregation: Compile a comprehensive quality report using MultiQC.
    3. Adapter Trimming: Remove adapter sequences and trim reads based on quality using bbduk.
    4. Post-trimming Quality Evaluation: Reassess the quality of reads post-trimming with FastQC.
    5. Alignment: Align the processed reads to a reference genome using the STAR aligner.
    6. Quantification: Count the number of reads mapping to genes using featureCounts.
Tools
    * Snakemake
    * FastQC
    * MultiQC
    * bbduk (part of the BBTools)
    * STAR aligner
    * samtools
    * featureCounts (part of the subread)
