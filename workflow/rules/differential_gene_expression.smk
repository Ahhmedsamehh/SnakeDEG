# =====================================================================
#  Stage 5 — differential gene expression with DESeq2  (the "DEG" rule)
# ---------------------------------------------------------------------
#  This rule uses `script:` instead of `shell:`. Snakemake injects an R
#  object called `snakemake` into the script, exposing snakemake@input,
#  @output, @params, @threads, @log — so we never parse argv by hand.
# =====================================================================

rule deseq2:
    input:
        # The raw count matrix from featureCounts ...
        counts=rules.featurecounts.output.counts,
        # ... and the sample sheet, which carries the condition labels.
        samples=config["samples"],
    output:
        # Full results for every gene (log2FC, p, padj).
        results="results/deg/deseq2_results.tsv",
        # Just the genes passing the padj + log2FC cutoffs.
        significant="results/deg/deseq2_significant.tsv",
        # DESeq2-normalised counts (useful for plotting / sharing).
        normalized="results/deg/normalized_counts.tsv",
        # Diagnostic + result figures.
        pca="results/deg/pca_plot.pdf",
        ma="results/deg/ma_plot.pdf",
        volcano="results/deg/volcano_plot.pdf",
        heatmap="results/deg/sample_distance_heatmap.pdf",
    params:
        # Everything the R script needs to know, pulled from config.
        condition_col=config["deseq2"]["condition_column"],
        reference_level=config["deseq2"]["reference_level"],
        treatment_level=config["deseq2"]["treatment_level"],
        padj=config["deseq2"]["padj_threshold"],
        lfc=config["deseq2"]["lfc_threshold"],
        # Sample names IN THE SAME ORDER featureCounts received the BAMs,
        # so we can rename the matrix's file-path columns reliably.
        sample_order=SAMPLES,
    log:
        "logs/deseq2.log",
    threads: 2
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2.R"
