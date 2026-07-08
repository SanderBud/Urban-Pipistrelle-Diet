
# Script to quality control raw seqs, clean/trim them, and assign taxonomy using a classifier (trained in 01_buildCOI_BOLD_NCBI_database.smk)
# For project Urban Bat Diet
# NIOO-KNAW
# By Judith E. Risse and Sander W. J. Buddendorf, 2025-2026


import os
import pandas as pd

base_dir = R"/path/to/server/data/result_X204SC24124169-Z01-F002/01.RawData"
sample_info_file = os.path.join(base_dir, "SampleSeq_info.csv")

# Take sample names from csv
df = pd.read_csv(sample_info_file, header=None)
samples = df.iloc[:, 0].astype(str).tolist()

output_dir_fastqc = os.path.join(base_dir, "fastqc_reports")
output_dir_fastp = R"/path/to/server/fastp"

container_qiime = "resources/amplicon_2024.10.sif"
container_fastqc = "resources/fastqc_0.12.1.sif"
container_multiqc = "resources/multiqc_1.26.sif"
container_fastp = "resources/fastp_0.24.0--heae3180_1.sif"

rule all:
    input:
        os.path.join(output_dir_fastqc, "multiqc_report.html"),
        expand(f"{output_dir_fastp}/{{sample}}/{{sample}}.clip_1.fastq.gz", sample=samples),
        expand(f"{output_dir_fastp}/{{sample}}/{{sample}}.clip_2.fastq.gz", sample=samples),
        "paired-end-demux.qza",
        "paired-end-trimmed.qza",
        "paired-end-trimmed.qzv",
        "fastq-stats.qzv",
        "dada2.stats.qzv",
        "dada2.repseqs.qzv",
        "dada2.table.qzv",
        "taxonomy_16S_SKLEARN_tempvar.qzv",
        "taxonomy_16S_SKLEARN_tempvar_barplot.qzv",
        "exported-feature-table/feature-table.tsv",
        "exported-taxonomy/taxonomy.tsv",
        "annotated_repseqs_16S_SKLEARN_tempvar.qzv",
        "taxonomy_16S_filteredMitoChloro_SKLEARN_tempvar_barplot.qzv",
        "taxonomy_16S_filteredMitoChloroEuk_SKLEARN_tempvar_barplot.qzv",
        "annotated_repseqs_16S_filteredMitoChloro_SKLEARN_tempvar.qzv",
        "annotated_repseqs_16S_filteredMitoChloroEuk_SKLEARN_tempvar.qzv"

""" This rule was needed because primers did not start at the first position 
    Trims the first 8 bases from forward and reverse reads (to remove primer offsets)
    Produces cleaned, paired FASTQ files """
rule fastp:
    input:
        raw1 = lambda wildcards: os.path.join(base_dir, wildcards.sample, f"{wildcards.sample}.raw_1.fastq.gz"),
        raw2 = lambda wildcards: os.path.join(base_dir, wildcards.sample, f"{wildcards.sample}.raw_2.fastq.gz"),

    output:
        clip1 = f"{output_dir_fastp}/{{sample}}/{{sample}}.clip_1.fastq.gz",
        clip2 = f"{output_dir_fastp}/{{sample}}/{{sample}}.clip_2.fastq.gz"
    params:
        clip="--trim_front1 8 --trim_front2 8" 
    singularity:
        container_fastp
    threads:
        20
    shell:
        """
            apptainer exec -B $PWD/cache:/path/to/server/fastp_0.24.0--heae3180_1.sif fastp -w {threads} --in1 {input.raw1} --in2 {input.raw2} --out1 {output.clip1} --out2 {output.clip2} {params.clip}
        """

""" Runs FastQC on raw (untrimmed) reads to generate per-sample quality reports """
rule fastqc:
    input:
        raw1 = lambda wildcards: os.path.join(base_dir, wildcards.sample, f"{wildcards.sample}.raw_1.fastq.gz"),
        raw2 = lambda wildcards: os.path.join(base_dir, wildcards.sample, f"{wildcards.sample}.raw_2.fastq.gz")
    output:
        report1 = os.path.join(output_dir_fastqc, "{sample}.raw_1_fastqc.html"),
        report2 = os.path.join(output_dir_fastqc, "{sample}.raw_2_fastqc.html")
    container: container_fastqc
    log:
        "logs/fastqc/{sample}.log"
    shell:
        """
        mkdir -p {output_dir_fastqc};
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 {container_fastqc} \
        fastqc -o {output_dir_fastqc} {input.raw1} {input.raw2} \
        2> {log}
        """

""" Aggregates all FastQC reports into a single MultiQC report """
rule multiqc:
    input:
        fastqc_reports = expand(os.path.join(output_dir_fastqc, "{sample}.raw_1_fastqc.html"), sample=samples) +
                         expand(os.path.join(output_dir_fastqc, "{sample}.raw_2_fastqc.html"), sample=samples)
    output:
        multiqc_report = os.path.join(output_dir_fastqc, "multiqc_report.html")
    container: container_multiqc
    log:
        "logs/multiqc.log"
    shell:
        """
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 {container_multiqc} \
        multiqc {output_dir_fastqc} -o {output_dir_fastqc} \
        2> {log}
        """

 """ Imports the fastp-trimmed reads into QIIME2 as a paired-end artifact (`.qza`) using a manifest file """
rule qiimeImport:
    input:
        "manifest.txt"
    output:
        "paired-end-demux.qza"
    container: container_qiime
    log:
        "logs/qiime_import.log"
    shell:
        """
            apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif \
            qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {input} \
            --output-path {output} \
            --input-format PairedEndFastqManifestPhred33V2 \
            2> {log}
        """

""" Trims remaining primer sequences using QIIME2's cutadapt plugin
    Discards reads where primers weren't found """
rule qiimeCutadapt:
    input:
        "paired-end-demux.qza",
    output:
        "paired-end-trimmed.qza",
    params:
        f="^GGDACWGGWTGAACWGTWTAYCCHCC...CCWGTWYTAGCHGGDGCWATYAC",
        r="^GTRATWGCHCCDGCTARWACWGG...GGDGGRTAWACWGTTCAWCCWGTHCC",
        extra="--p-discard-untrimmed" 
    log:
        "logs/qiime_trimmed.log"
    container: container_qiime
    threads:
        40 
    shell:
        "apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif qiime cutadapt trim-paired --i-demultiplexed-sequences {input} --p-cores {threads} --p-adapter-f {params.f} --p-adapter-r {params.r} --o-trimmed-sequences {output} {params.extra} --verbose &> {log}"

""" Generates a summary visualization (`.qzv`) of sequence counts per sample after primer trimming """
rule qiimeTrimmedMetrics:
    input:
        "paired-end-trimmed.qza"
    output:
        "paired-end-trimmed.qzv"
    log:
        "logs/qiime_tMetrics.log"
    container: "resources/amplicon_2024.10.sif"
    shell:
        """
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif \
        qiime demux summarize \
        --i-data {input} \
        --o-visualization {output} \
        2> {log} 
        """

""" Calculates quality statistics of trimmed reads (length, quality distribution) using vsearch """
rule qiimeVsearchStats:
    input:
        "paired-end-trimmed.qza"
    output:
        "fastq-stats.qzv"
    log:
        "logs/qiime_fastqstats.log"
    container: "resources/amplicon_2024.10.sif"
    threads:
        10
    shell:
        """  
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif \
        qiime vsearch fastq-stats \
        --i-sequences {input} \
        --p-threads {threads} \
        --o-visualization {output}  \
        --verbose &> {log}
        """

""" Performs denoising, chimera removal, and dereplication using DADA2
    Outputs the representative ASV sequences, the feature table, and denoising stats """
rule qiimeDADA2:
    input:
        "paired-end-trimmed.qza"
    output:
        repseqs="dada2.repseqs.qza",
        table="dada2.table.qza",
        stats="dada2.stats.qza"
    params:
        trimLeft="0",
        truncLeft="150",
        trimRight="0",
        truncRight="150",
        maxEEleft=2,
        maxEEright=2
    log:
        "logs/qiime_dada2.log"
    threads:
        30
    container: "resources/amplicon_2024.10.sif"
    shell:
        """
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif \
        qiime dada2 denoise-paired \
        --p-n-threads {threads} \
        --i-demultiplexed-seqs {input} \
        --o-representative-sequences {output.repseqs} \
        --o-table {output.table} \
        --o-denoising-stats {output.stats} \
        --p-trim-left-f {params.trimLeft} \
        --p-trunc-len-f {params.truncLeft} \
        --p-trim-left-r {params.trimRight} \
        --p-trunc-len-r {params.truncRight} \
        --p-max-ee-f {params.maxEEleft} \
        --p-max-ee-r {params.maxEEright} \
        --verbose &> {log}
        """

""" Visualizes DADA2’s denoising statistics (e.g., how many reads were retained per step) """
rule qiimeDADA2Metrics:
    input:
        "dada2.stats.qza"
    output:
        "dada2.stats.qzv"
    log:
        "logs/qiime_dada2Metrics.log"
    container: "resources/amplicon_2024.10.sif"
    shell:
        """
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif \
        qiime metadata tabulate \
        --m-input-file {input} \
        --o-visualization {output}
        """

""" Summarizes the ASV feature table — i.e., shows how many reads per sample and per ASV """
rule qiimeDADA2Summarize:
    input:
        dada2="dada2.table.qza",
        sampleData="manifest.txt"
    output:
        "dada2.table.qzv"
    log:
        "logs/qiime_dada2Summarize.log"
    container: "resources/amplicon_2024.10.sif"
    shell:
        """
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif \
        qiime feature-table summarize \
        --i-table {input.dada2} \
        --o-visualization {output} \
        --m-sample-metadata-file {input.sampleData} \
        2>{log}
        """

""" Tabulates and visualizes the ASV sequences (fasta-style) in an interactive `.qzv` file """
rule qiimeDADA2Tabulate:
    input:
        dada2="dada2.repseqs.qza",
    output:
        "dada2.repseqs.qzv"
    log:
        "logs/qiime_dada2Tabulate.log"
    container: "resources/amplicon_2024.10.sif"
    shell:
        """
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif \
        qiime feature-table tabulate-seqs \
        --i-data {input} \
        --o-visualization {output} \
        2>{log}
        """

""" Assigns taxonomy to each ASV using a pre-trained Naive Bayes classifier (trained on BOLD and NCBI)
    Outputs a `.qza` taxonomy file with confidence scores """
rule qiimeFeatureClassifier:
    input:
        repseqs="dada2.repseqs.qza",
    output:
        tax_temp="taxonomy_16S_SKLEARN_tempvar.qza"
    params:
        classifier="fitClassifier_boldNCBI_COI/classifier.qza"
    log:
        "logs/qiime_FeatureClassifier.log"
    threads:
        20
    container: "resources/amplicon_2024.10.sif"
    shell:
        "apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif qiime feature-classifier classify-sklearn --i-classifier {params.classifier} --i-reads {input.repseqs} --o-classification {output.tax_temp} --p-n-jobs {threads} 2>{log}"

""" Tabulates the taxonomy classification output for visualization. """
rule qiimeFeatureClassifierTabulate:
        input:
            "taxonomy_16S_SKLEARN_tempvar.qza",
        output:
            "taxonomy_16S_SKLEARN_tempvar.qzv"
        log:
            "logs/qiime_FeatureClassifierTabulate.log"
        container: "resources/amplicon_2024.10.sif"
        shell:
           "apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif qiime metadata tabulate   --m-input-file {input} --o-visualization {output} 2>{log}"

""" Tabulates and visualizes both taxonomy and ASV sequences. """
rule qiimeRepSeqsTabulate:
        input:
            tax="taxonomy_16S_SKLEARN_tempvar.qza",
            repseqs="dada2.repseqs.qza"
        output:
            "annotated_repseqs_16S_SKLEARN_tempvar.qzv"
        log:
            "logs/qiime_FeatureClassifierTabulate.log"
        container: "resources/amplicon_2024.10.sif"
        shell:
           "apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif qiime metadata tabulate   --m-input-file {input.tax} --m-input-file {input.repseqs}  --o-visualization {output} 2>{log}"

""" Creates a barplot of taxa abundances for each sample. """
rule qiimeTaxaBarplot:
        input:
            table="dada2.table.qza",
            taxonomy="taxonomy_16S_SKLEARN_tempvar.qza",
            sampleData="manifest.txt"
        output:
            "taxonomy_16S_SKLEARN_tempvar_barplot.qzv"
        log:
            "logs/qiime_TaxaBarplot.log"
        container: "resources/amplicon_2024.10.sif"
        shell:
            "apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif qiime taxa barplot --i-table {input.table} --i-taxonomy {input.taxonomy} --m-metadata-file {input.sampleData} --o-visualization {output}"

""" Exports the taxonomy """
rule export_taxonomy_tsv:
    input:
        tax="taxonomy_16S_SKLEARN_tempvar.qza"
    output:
        tsv="exported-taxonomy/taxonomy.tsv"
    shell:
        """
        mkdir -p exported-taxonomy
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif \
            qiime tools export --input-path {input.tax} --output-path exported-taxonomy
        """

""" Exports the feature table """
rule export_feature_table_tsv:
    input:
        table="dada2.table.qza"
    output:
        tsv="exported-feature-table/feature-table.tsv"
    shell:
        """
        mkdir -p exported-feature-table
        apptainer exec -B $PWD/cache:/path/to/server/qiime2 resources/amplicon_2024.10.sif \
            qiime tools export --input-path {input.table} --output-path exported-feature-table
        
        # Convert biom to tsv, input biom is inside exported-feature-table/feature-table.biom
        biom convert -i exported-feature-table/feature-table.biom -o {output.tsv} --to-tsv --header-key taxonomy
        """
