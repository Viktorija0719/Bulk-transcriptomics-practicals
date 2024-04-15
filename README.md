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

### Initial Directory Structure

The project directory is organized as follows:

- `data/`: This directory contains all data-related directories.
  - `data_ref_annot/`: Inside this directory, you'll find reference and annotation files necessary for the workflow.
    - `adapters.fa`: Adapter sequences for trimming.
    - `chr19_20Mb.bed`: BED file for the genomic region of interest.
    - `chr19_20Mb.fa`: Reference genome sequence for the specified genomic region.
    - `chr19_20Mb.gtf`: Gene annotation file in GTF format.
  - `RNA-seq-sample-data/`: This directory contains all the RNA-seq sample data files.


### Workflow-Generated Directories

The workflow will create and use the following directories:

- `data/processed`: Stores trimmed sequencing reads.
- `quality/raw_fastqc` and `quality/proc_fastqc`: Contain FastQC reports for raw and trimmed reads, respectively.
- `quality/raw_multiqc` and `quality/proc_multiqc`: Contain MultiQC reports aggregating FastQC results before and after trimming.
- `mapping/mapping_1`: Stores STAR aligner output files.
- `mapping/sort_1`: Contains sorted BAM files.
- `count/counts`: Stores gene quantification results from featureCounts.


### Running the Workflow


![image](https://github.com/Viktorija0719/Bulk-transcriptomics-practicals/assets/150614034/eeac6093-e0fe-4789-8852-1093cff44369)

![image](https://github.com/Viktorija0719/Bulk-transcriptomics-practicals/assets/150614034/eea18a50-a1cf-4223-b1a9-085ce4b79f05)


