# =====================================================================
#  Stage 4 — gene-level read counting with featureCounts
# ---------------------------------------------------------------------
#  Unlike the per-sample steps above, this rule has NO {sample} wildcard:
#  it takes ALL sorted BAMs at once and emits a single count matrix
#  (genes x samples) — exactly the input DESeq2 expects.
# =====================================================================

rule featurecounts:
    input:
        # expand() builds the full BAM list. The ORDER here defines the
        # column order in the matrix, and we reuse SAMPLES (same order)
        # in the DESeq2 script to label those columns correctly.
        bams=expand("results/aligned/{sample}.sorted.bam", sample=SAMPLES),
        gtf=config["reference"]["gtf"],
    output:
        # featureCounts always writes "<counts>.summary" alongside the
        # matrix; we declare it so MultiQC can pick it up.
        counts="results/counts/all_samples_counts.tsv",
        summary="results/counts/all_samples_counts.tsv.summary",
    params:
        strandedness=config["featurecounts"]["strandedness"],
    log:
        "logs/featurecounts.log",
    threads: config["threads"]
    conda:
        "../envs/subread.yaml"
    shell:
        # -p --countReadPairs : count fragments (read pairs), not reads.
        # -t exon -g gene_id  : sum exon counts up to the gene level.
        # -s {strandedness}   : 0/1/2 from config (see config.yaml).
        """
        featureCounts \
            -T {threads} \
            -p --countReadPairs \
            -t exon -g gene_id \
            -s {params.strandedness} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} \
            > {log} 2>&1
        """
