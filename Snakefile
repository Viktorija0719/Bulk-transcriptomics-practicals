configfile: "config.yaml"

import os
os.makedirs("data/processed", exist_ok=True)
os.makedirs("quality/raw_fastqc", exist_ok=True)
os.makedirs("quality/raw_multiqc", exist_ok=True)
os.makedirs("quality/proc_fastqc", exist_ok=True)
os.makedirs("data/indexed", exist_ok=True)
os.makedirs("mapping/mapping_1", exist_ok=True)
os.makedirs("mapping/index_1", exist_ok=True)
os.makedirs("mapping/sort_1", exist_ok=True)
os.makedirs("count/", exist_ok=True)

rule all:
    input:
        expand("quality/proc_multiqc/multiqc_trimmed_report.html", sample=config['samples']),
        expand("quality/raw_multiqc/multiqc_report.html", sample=config['samples']),
        expand("mapping/mapping_1/{sample}.Aligned.sortedByCoord.out.bam", sample=config['samples']),
        expand("mapping/index_1/{sample}.sorted.bam.bai", sample=config['samples']),
        expand("count/counts/{sample}_counts.txt", sample=config['samples'])


rule fastqc_raw:
    input:
        r1="data/RNA-seq-sample-data/{sample}_R1_001.fastq.gz",
        r2="data/RNA-seq-sample-data/{sample}_R2_001.fastq.gz"
    output:
        r1_html="quality/raw_fastqc/{sample}_R1_001_fastqc.html",
        r1_zip="quality/raw_fastqc/{sample}_R1_001_fastqc.zip",
        r2_html="quality/raw_fastqc/{sample}_R2_001_fastqc.html",
        r2_zip="quality/raw_fastqc/{sample}_R2_001_fastqc.zip"
    shell:
        "fastqc {input.r1} {input.r2} --outdir quality/raw_fastqc"

rule multiqc_raw:
    input:
        lambda wildcards: expand("quality/raw_fastqc/{sample}_{read}_001_fastqc.zip", sample=config['samples'], read=['R1', 'R2'])
    output:
        "quality/raw_multiqc/multiqc_report.html"
    shell:
        "multiqc quality/raw_fastqc/ --filename {output} -f"


rule bbduk_trim:
    input:
        r1="data/RNA-seq-sample-data/{sample}_R1_001.fastq.gz",
        r2="data/RNA-seq-sample-data/{sample}_R2_001.fastq.gz",
        adapters=config['directories']['adapters']
    output:
        trimmed_r1="data/processed/{sample}_R1_001_trimmed.fastq.gz",
        trimmed_r2="data/processed/{sample}_R2_001_trimmed.fastq.gz"
    shell:
        """
        bbduk.sh \
            in1={input.r1} \
            in2={input.r2} \
            out1={output.trimmed_r1} \
            out2={output.trimmed_r2} \
            ref={input.adapters} \
            ktrim=r \
            k=23 \
            mink=11 \
            hdist=1 \
            tpe=t \
            tbo=t \
            qtrim=r \
            trimq=10
        """

rule fastqc_processed:
    input:
        r1="data/processed/{sample}_R1_001_trimmed.fastq.gz",
        r2="data/processed/{sample}_R2_001_trimmed.fastq.gz"
    output:
        r1_html="quality/proc_fastqc/{sample}_R1_001_trimmed_fastqc.html",
        r1_zip="quality/proc_fastqc/{sample}_R1_001_trimmed_fastqc.zip",
        r2_html="quality/proc_fastqc/{sample}_R2_001_trimmed_fastqc.html",
        r2_zip="quality/proc_fastqc/{sample}_R2_001_trimmed_fastqc.zip"
    shell:
        "fastqc {input.r1} {input.r2} --outdir quality/proc_fastqc"



rule multiqc_processed:
    input:
        lambda wildcards: expand("quality/proc_fastqc/{sample}_{read}_001_trimmed_fastqc.zip", sample=config['samples'], read=['R1', 'R2'])
    output:
        "quality/proc_multiqc/multiqc_trimmed_report.html"
    shell:
        "multiqc quality/proc_fastqc/ --filename multiqc_trimmed_report.html -o quality/proc_multiqc/"


rule star_genome_index:
    input:
        fasta="data/data_ref_annot/chr19_20Mb.fa",
        gtf="data/data_ref_annot/chr19_20Mb.gtf"
    output:
        directory("data/indexed/exonInfo.tab")
    params:
        genomeSAindexNbases=11
    threads: 16
    shell:
        """
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --genomeSAindexNbases {params.genomeSAindexNbases}
        """



rule star_align:
    input:
        r1="data/processed/{sample}_R1_001_trimmed.fastq.gz",
        r2="data/processed/{sample}_R2_001_trimmed.fastq.gz",
        index_dir=("data/indexed/exonInfo.tab")
    output:
        bam="mapping/mapping_1/{sample}.Aligned.sortedByCoord.out.bam"
    threads: config['star']['runThreadN']
    shell:
        """
        STAR --genomeDir {input.index_dir} \
             --readFilesIn {input.r1} {input.r2} \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix mapping/mapping_1/{wildcards.sample}. \
             --runThreadN {threads}
        """

# Sort BAM files using samtools
rule samtools_sort:
    input:
        "mapping/mapping_1/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        sorted_bam="mapping/index_1/{sample}.sorted.bam"
    shell:
        "samtools sort -o {output.sorted_bam} {input}"

# Index sorted BAM files
rule samtools_index:
    input:
        "mapping/index_1/{sample}.sorted.bam"
    output:
        "mapping/index_1/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"


def get_strand_specificity(sample):
    if "_R1_" in sample:
        return "1"  # Forward strand
    elif "_R2_" in sample:
        return "2"  # Reverse strand
    else:
        return "0"  # Unstranded


rule featureCounts:
    input:
        bam="mapping/index_1/{sample}.sorted.bam",  # Adjusted to match the sorted BAM file path
        annotation=config['featureCounts']['annotation']
    output:
        counts="count/counts/{sample}_counts.txt"
    params:
        strand_specificity=lambda wildcards: get_strand_specificity(wildcards.sample),
    shell:
        """
        featureCounts -a {input.annotation} -o {output.counts} -s {params.strand_specificity} -p -t exon -g gene_id {input.bam}
        """

