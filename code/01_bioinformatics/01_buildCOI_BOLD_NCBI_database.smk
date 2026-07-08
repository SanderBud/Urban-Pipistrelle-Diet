
# Script to clean BOLD COI data, download and clean NCBI COI data, and train a classifier.
# For project Urban Bat Diet
# NIOO-KNAW
# By Judith E. Risse, 2025-2026


SAMPLES, = glob_wildcards("/path/to/server/boldNCBI_COI_Classifier/boldCustom.{sample}.seqNtaxa.csv")

QIIME_CONTAINER = "/path/to/server/software/singularity/amplicon_2024.10.sif"
SEQKIT_CONTAINER = "/path/to/server/software/singularity/seqkit_2.2.0--h9ee0642_0.sif"
MAFFT_CONTAINER = "/path/to/server/software/singularity/mafft_7.490--hec16e2b_1.sif"

rule all:
    input:
        "COI_primers.fasta.map_toSmall",
        "COI_primers.fasta.map_toGiant",
        "tmpdir_boldNCBIFullSeqs/taxonomy.tsv",
        "fitClassifier_boldNCBI_COI/classifier.qza",
        "fitClassifier_boldNCBI_COI/evaluation.qzv",
        "fitClassifier_boldNCBI_COI/observed_taxonomy.qza",

rule concatenateInput:
    input:
        expand("/path/to/server/boldCustom.{sample}.seqNtaxa.csv", sample=SAMPLES)
    output:
        "bold_allrawSeqs.fasta",
    shell:
        """cat /path/to/server/boldNCBI_COI_Classifier/*.seqNtaxa.csv | grep -v 'sequenceID,taxon,nucleotides' | sort | uniq | awk -F ',' '{{OFS="\\t"}};{{print $1";tax="$2, $3}}' | sed 's/^/>/g' | tr '\\t' '\\n' > {output}"""

rule removeNonIUPAC:
    input:
        "bold_allrawSeqs.fasta",
    output:
        "cleaned_bold_allrawSeqs.fasta",   
    shell:
        """
            ## getting the line numbers of the bad sequences
            grep -v '^>' {input} | grep [^THADGRC\.SMBNWVKY-] -n > badseqs.txt

            ## reformatting these line numbers to delete
            for val in $(cut -f 1 badseqs.txt -d ':'); do
            myvar1=$(($val*2))
            myvar2=$(($myvar1-1))
            echo $myvar2"d"
            echo $myvar1"d"
            done | tr '\\n' ';' | sed 's/.$//' > droplines.txt
            
            ## deleting these line numbers
            myval=$(cat droplines.txt)
            sed -e $myval {input} > {output}
        """

rule makeQiimeInput:
    input:
         "cleaned_bold_allrawSeqs.fasta",
    output:
        seqFile="bold_rawSeqs_forQiime.fasta",
        taxFile="bold_rawTaxa_forQiime.tsv",
    shell:
        """
            ## First make the fasta file
            sed 's/;tax=k.*//' {input} | paste - - | awk '{{gsub("-","",$2);print($1"\\n"$2)}}'> {output.seqFile}
            
            ## Next make the taxonomy file
            grep '^>' {input} | sed 's/^>//' |awk '{{for(i=1;i<2;i++)sub(";","\\t")}}1' > {output.taxFile}
        """

rule importQiimeSeq:
    input:
        seqFile="bold_rawSeqs_forQiime.fasta",
    output:
        "bold_rawSeqs.qza",
    container:
        QIIME_CONTAINER,
    shell:
        """
            qiime tools import --input-path {input.seqFile} --output-path {output} --type 'FeatureData[Sequence]'
        """

rule importQiimeTax:
    input:
        taxFile="bold_rawTaxa_forQiime.tsv",
    output:
        "bold_rawTaxa.qza"
    container:
        QIIME_CONTAINER,
    shell:
        """qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path {input.taxFile} --output-path {output} """

rule downloadNCBI:
    output:
        NCBI="NCBIdata_notBOLD_sequences.qza",
        NCBItax="NCBIdata_notBOLD_taxonomy.qza",
    container:
        QIIME_CONTAINER,
    threads:
        10
    shell:
        """
            qiime rescript get-ncbi-data --p-query '(cytochrome c oxidase subunit 1[Title] OR cytochrome c oxidase subunit I[Title] OR cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title] OR COX1[Title] OR CO1[Title] OR COI[Title]) AND mitochondrion[Filter] NOT environmental sample[Title] NOT environmental samples[Title] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title] NOT "BARCODE"[KYWD])' --o-sequences {output.NCBI} --o-taxonomy {output.NCBItax} --p-n-jobs {threads}
        """


rule qiimeMergeNCBI:
    input:
        NCBI="NCBIdata_notBOLD_sequences.qza",
        NCBItax="NCBIdata_notBOLD_taxonomy.qza",
        BOLD="bold_rawSeqs.qza",
        BOLDtax="bold_rawTaxa.qza",
    output:
        tax="boldNCBI_rawTaxa.qza",
        seqs="boldNCBI_rawSeqs.qza",
    container:
        QIIME_CONTAINER,
    shell:
        """ 
            qiime rescript merge-taxa --i-data {input.BOLDtax}  {input.NCBItax} --o-merged-data {output.tax};
            qiime feature-table merge-seqs --i-data {input.BOLD} {input.NCBI} --o-merged-data {output.seqs}
        """

rule qiimeDegap1:
    input:
       "bold_rawSeqs.qza",
    output:
        "bold_degapSeqs_min1.qza", 
    container:
        QIIME_CONTAINER,
    shell:
        "qiime rescript degap-seqs --i-aligned-sequences {input} --p-min-length 1 --o-degapped-sequences {output}"

rule exportSeqs:
    input:
        "boldNCBI_rawSeqs.qza",
    output:
        "degap_tmpdir/dna-sequences.fasta",
    container:
        QIIME_CONTAINER,
    shell:
        """
            qiime tools export --input-path {input} --output-path degap_tmpdir
        """


rule trimSeqs:
    input:
        "degap_tmpdir/dna-sequences.fasta",
    output:
        "boldNCBI_outerNtrimmed.fasta", 
    shell:
        """
            cat {input} | paste - - | awk -v OFS="\\t" '{{gsub(/^[N]+/,"",$2);print}}' | rev | sed 's/^N*//' | rev | tr '\\t' '\\n' > {output}
        """

rule importQiimeSeq2:
    input:
        seqFile="boldNCBI_outerNtrimmed.fasta",
    output:
        "boldNCBI_outerNtrimmed.qza",
    container:
        QIIME_CONTAINER,
    shell:
        """
            qiime tools import --input-path {input.seqFile} --output-path {output} --type 'FeatureData[Sequence]'"""

        

rule qiimeCullSeqs:
    input:
       "boldNCBI_outerNtrimmed.qza",
    output:
        "boldNCBI_ambi_hpoly_filtd_seqs.qza", 
    threads:
        20
    container:
        QIIME_CONTAINER,
    shell:
        """qiime rescript cull-seqs --i-sequences {input} --p-n-jobs {threads} --p-num-degenerates 5 --p-homopolymer-length 12 --o-clean-sequences {output} """

rule qiimeFilterLength:
    input:
       "boldNCBI_ambi_hpoly_filtd_seqs.qza",
    output:
        kept="boldNCBI_ambi_hpoly_length_filtd_seqs.qza",
        disc="boldNCBI_ambi_hpoly_length_discarded_seqs.qza", 
    container:
        QIIME_CONTAINER,
    shell:
        """qiime rescript filter-seqs-length --i-sequences {input} --p-global-min 250 --p-global-max 1600 --o-filtered-seqs {output.kept} --o-discarded-seqs {output.disc}"""

rule qiimeDerep1:
    input:
        seqs="boldNCBI_ambi_hpoly_length_filtd_seqs.qza",
        taxa="boldNCBI_rawTaxa.qza",
    output:
        seqs="boldNCBI_derep1_seqs.qza",
        taxa="boldNCBI_derep1_taxa.qza"
    container:
        QIIME_CONTAINER
    shell:
        """
        qiime rescript dereplicate \
        --i-sequences {input.seqs} \
        --i-taxa {input.taxa} \
        --p-mode 'super' \
        --p-derep-prefix \
        --o-dereplicated-sequences {output.seqs} \
        --o-dereplicated-taxa {output.taxa}
        """


rule exportSeqsDerep:
    input:
        "boldNCBI_derep1_seqs.qza",
    output:
        "tmpdir_boldNCBIFullSeqs/dna-sequences.fasta",
    params:
        "tmpdir_boldNCBIFullSeqs"
    container:
        QIIME_CONTAINER,
    shell:
        """
            qiime tools export --input-path {input} --output-path {params}
        """

rule exportTaxaDerep:
    input:
        "boldNCBI_derep1_taxa.qza",
    output:
        "tmpdir_boldNCBIFullSeqs/taxonomy.tsv",
    params:
        "tmpdir_boldNCBIFullSeqs"
    container:
        QIIME_CONTAINER,
    shell:
        """
            qiime tools export --input-path {input} --output-path {params}
        """

rule seqkitTrimLength:
    input:
        "tmpdir_boldNCBIFullSeqs/dna-sequences.fasta",
    output:
        "boldNCBIFull_lengthFiltd_seqs.fasta",
    container:
        SEQKIT_CONTAINER
    shell:
        "seqkit seq --min-len 660 --max-len 1000 -w 0 {input} > {output}"

rule seqkitPrepMafftSmall: 
    input:
        tax="tmpdir_boldNCBIFullSeqs/taxonomy.tsv",
        seq="boldNCBIFull_lengthFiltd_seqs.fasta",
    output:
        "boldNCBI_mafftRefs_subseqs_seqs.fasta"
    container:
        SEQKIT_CONTAINER
    shell:
        """
            grep 'p__Arthropoda;c__Insecta' {input.tax} | grep -v 's__$' | grep -v 's__.*sp.$' | grep -v '-' | cut -f 1 > insect_species.list;
            seqkit grep --pattern-file insect_species.list {input.seq} > boldFull_lengthNtaxaFiltd_seqs.fasta;
            seqkit fx2tab -n -i -a boldFull_lengthNtaxaFiltd_seqs.fasta | awk '{{if ($2 == "ACGT") print $1}}' > nonAmbig_featureIDs.txt;
            seqkit grep --pattern-file nonAmbig_featureIDs.txt boldFull_lengthNtaxaFiltd_seqs.fasta -w 0 > boldFull_lengthNtaxaNambigFiltd_seqs.fasta;
            seqkit sample --rand-seed 103 --number 3000 --two-pass boldFull_lengthNtaxaNambigFiltd_seqs.fasta -w 0 > {output};
        """

rule mafftSmall:
    input:
        seqs="boldNCBI_mafftRefs_subseqs_seqs.fasta",
        primers="COI_primers.fasta",
    output:
        primer="ref_primer.msa",
        pmap="COI_primers.fasta.map_toSmall",
    container:
        MAFFT_CONTAINER,
    threads:
        40
    shell:
        """
            mafft --auto --thread -1 {input.seqs} > reference_MSA;
            mafft --multipair --addfragments {input.primers} --keeplength --thread {threads} --mapout --reorder reference_MSA > {output.primer};
            mv {input.primers}.map {output.pmap}
        """

rule extractRegion:
    input:
        "{prefix}.msa",
    output:
        "{prefix}_coordinateTrimmed.fasta",
    shell:
        " ~/software/python3.6_venv/bin/python ../Sander/extract_alignment_region.py -i {input} -o {output} -s 1660  -e 2008 "

rule seqkitPrepMafftLarge:
    input: 
        refp="ref_primer_coordinateTrimmed.fasta",
        seq="tmpdir_boldNCBIFullSeqs/dna-sequences.fasta",
        maffts="boldNCBI_mafftRefs_subseqs_seqs.fasta",
    output:
        "seqsForGiantAlignment.fasta",
    container:
        SEQKIT_CONTAINER
    shell:
        """ 
            seqkit seq -g -w0 {input.refp} | paste - - | awk '{{print length ($2)}}' | sort | uniq -c;
            grep '^>' {input.maffts} | sed 's/^>//' > droplist.txt;
            seqkit grep -v --pattern-file droplist.txt {input.seq} > {output};
        """

rule mafftBig:
    input:
        seqs="seqsForGiantAlignment.fasta",
        align="ref_primer.msa",
    output:
        "boldNCBI_derep1_COI_tmp.msa",
    container:
        MAFFT_CONTAINER
    threads:
        40
    shell:
        "mafft --auto --addfull {input.seqs} --keeplength --thread {threads} {input.align} > {output}"

rule checkBig:
    input:
        msa="boldNCBI_derep1_COI_tmp.msa",
        primers="COI_primers.fasta",
    output:
        "COI_primers.fasta.map_toGiant",
    threads:
        40
    shell:
        """
            singularity exec /path/to/server/software/singularity/seqkit_2.2.0--h9ee0642_0.sif seqkit head -n 2000 {input.msa} > subset_giant_alignment
            singularity exec /path/to/server/software/singularity/mafft_7.490--hec16e2b_1.sif mafft --multipair --addfragments {input.primers} --keeplength --thread {threads} --mapout subset_giant_alignment > subset_giant_alignment_secondaryPrimerAlign
            mv COI_primers.fasta.map {output}
        """

rule removePrimers:
    input:
        "boldNCBI_derep1_COI_tmp_coordinateTrimmed.fasta",
    output:
        fasta="boldNCBI_derep1_COI.fasta",
        freq="COI_seq_lengths.freq.table",
    container: 
        SEQKIT_CONTAINER
    shell:
        """
            grep PRIMER {input} | sed 's/>//' | sed 's/[ \\t]*$//' > dropPrimers.txt
            seqkit grep -v -f dropPrimers.txt {input} | seqkit seq --upper-case -w 0 > {output.fasta}
            grep -v '^>' {output.fasta} | sed 's/-//g' | awk '{{print length}}' | sort | uniq -c > {output.freq}
        """

rule qiimeImportFinal:
    input:
        "boldNCBI_derep1_COI.fasta",
    output:
        "boldNCBI_derep1_COI_seq_aligned.qza",
    container:
        QIIME_CONTAINER,
    shell:
        "qiime tools import --input-path {input} --output-path {output} --type 'FeatureData[AlignedSequence]'"

rule qiimeDegapFinal:
    input:
        "boldNCBI_derep1_COI_seq_aligned.qza"
    output:
        "boldNCBI_derep1_COI_seqs_nogaps.qza",
    container:
        QIIME_CONTAINER,
    shell:
        "qiime rescript degap-seqs --i-aligned-sequences {input} --p-min-length 200 --o-degapped-sequences {output}"

rule qiimeDerepFinal:
    input:
        seqs="boldNCBI_derep1_COI_seqs_nogaps.qza", 
        tax="boldNCBI_rawTaxa.qza"
    output:
        tax="boldNCBI_COI_taxa.qza",
        seqs="boldNCBI_COI_seqs.qza",
    container:
        QIIME_CONTAINER,
    shell:
        "qiime rescript dereplicate --i-sequences {input.seqs} --i-taxa {input.tax} --p-mode 'super' --p-derep-prefix --o-dereplicated-sequences {output.seqs} --o-dereplicated-taxa {output.tax}"

rule qiimeTrainClassifier:
    input:
        tax="boldNCBI_COI_taxa.qza",
        seqs="boldNCBI_COI_seqs.qza",
    output:
        clas="fitClassifier_boldNCBI_COI/classifier.qza",
        eval="fitClassifier_boldNCBI_COI/evaluation.qzv",
        tax="fitClassifier_boldNCBI_COI/observed_taxonomy.qza",
    container:
        QIIME_CONTAINER,
    params:
        rpb=6000,
        outdir="fitClassifier_boldNCBI_COI",
    threads:
        10
    shell:
        "qiime rescript evaluate-fit-classifier --i-sequences {input.seqs} --i-taxonomy {input.tax} --p-reads-per-batch {params.rpb} --p-n-jobs {threads} --o-classifier {output.clas} --o-evaluation {output.eval} --o-observed-taxonomy {output.tax}"

