# =====================================================================
#  Stage 3 — alignment to the genome with HISAT2
# ---------------------------------------------------------------------
#  Two rules:
#    hisat2_index  : build the genome index ONCE (no {sample} wildcard).
#    hisat2_align  : map each sample, pipe straight to a sorted+indexed
#                    BAM (no giant intermediate .sam on disk).
# =====================================================================

rule hisat2_index:
    # Build the index from the reference FASTA named in config.yaml.
    input:
        fasta=config["reference"]["fasta"],
    output:
        # HISAT2 writes 8 files: <prefix>.1.ht2 ... <prefix>.8.ht2.
        # multiext() is shorthand for "same prefix, these suffixes".
        # Listing them as outputs lets Snakemake know the index is built.
        multiext(
            config["hisat2_index_prefix"],
            ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
            ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2",
        ),
    params:
        # The shell needs the bare prefix, not the file list.
        prefix=config["hisat2_index_prefix"],
    log:
        "logs/hisat2_index.log",
    threads: config["threads"]
    conda:
        "../envs/align.yaml"
    # NOTE: indexing a full human genome needs a LOT of RAM (~6 GB for
    # the default index; >160 GB if you ever add SNPs/transcripts) and
    # can take an hour+. It only runs once, then every sample reuses it.
    shell:
        """
        hisat2-build -p {threads} {input.fasta} {params.prefix} > {log} 2>&1
        """


rule hisat2_align:
    input:
        r1=rules.fastp.output.r1,
        r2=rules.fastp.output.r2,
        # Depend on the index files so this waits for hisat2_index.
        index=rules.hisat2_index.output,
    output:
        bam="results/aligned/{sample}.sorted.bam",
        bai="results/aligned/{sample}.sorted.bam.bai",
    params:
        prefix=config["hisat2_index_prefix"],
    log:
        # Two separate logs. The hisat2 one holds the alignment-rate
        # summary that MultiQC parses, so keep it on its own.
        hisat2="logs/hisat2/{sample}.log",
        sort="logs/samtools/{sample}.log",
    threads: config["threads"]
    conda:
        "../envs/align.yaml"
    shell:
        # --dta tunes splicing reporting for downstream assembly/counting.
        # The `|` streams SAM from hisat2 into samtools sort with no temp
        # SAM file; then we index the sorted BAM.
        """
        hisat2 -p {threads} --dta \
            -x {params.prefix} \
            -1 {input.r1} -2 {input.r2} 2> {log.hisat2} \
        | samtools sort -@ {threads} -o {output.bam} - 2> {log.sort}
        samtools index {output.bam}
        """
