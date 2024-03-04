# RNA-seq Analysis Workflow

## Project Description

This project contains a Snakemake workflow designed for comprehensive RNA sequencing (RNA-seq) data analysis. It includes key steps such as quality control, adapter trimming, read alignment, and quantification.

## Workflow Steps

1. **Quality Evaluation**: Assess the quality of raw sequencing reads using FastQC.
2. **Quality Report Aggregation**: Compile a comprehensive quality report using MultiQC.
3. **Adapter Trimming**: Remove adapter sequences and trim reads based on quality using bbduk.
4. **Post-trimming Quality Evaluation**: Reassess the quality of reads post-trimming with FastQC.
5. **Alignment**: Align the processed reads to a reference genome using the STAR aligner.
6. **Quantification**: Count the number of reads mapping to genes using featureCounts.

## Tools

- **Snakemake**: A workflow management system that helps to create and manage complex bioinformatics workflows.
- **FastQC**: A quality control tool for high throughput sequence data.
- **MultiQC**: A tool to aggregate bioinformatics results across many samples into a single report.
- **bbduk (part of BBTools)**: A tool for trimming adapters and other contaminants from sequencing reads.
- **STAR aligner**: A fast RNA-seq read mapper.
- **samtools**: A suite of programs for interacting with high-throughput sequencing data.
- **featureCounts (part of subread)**: A software package used for counting reads to genomic features such as genes.

## Usage

Your data should be organized as follows for the workflow to function correctly:
data/
├── data_ref_annot/
│   ├── adapters.fa
│   ├── chr19_20Mb.bed
│   ├── chr19_20Mb.fa
│   └── chr19_20Mb.gtf
└── RNA-seq-sample-data/
    ├── Sample_1_R1.fastq.gz
    ├── Sample_1_R2.fastq.gz
    ...



```markdown
```
.
├── data/
│   ├── data_ref_annot/
│   │   ├── adapters.fa
│   │   ├── chr19_20Mb.bed
│   │   ├── chr19_20Mb.fa
│   │   └── chr19_20Mb.gtf
│   └── RNA-seq-sample-data/
│       ├── Sample_1_R1.fastq.gz
│       ├── Sample_1_R2.fastq.gz
│       ...
```
```
