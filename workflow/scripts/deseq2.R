# =====================================================================
#  deseq2.R  —  differential expression for the Snakemake "deseq2" rule
# ---------------------------------------------------------------------
#  Snakemake runs this with `script:`, which means it injects an S4
#  object called `snakemake`. We read everything from it:
#     snakemake@input[["counts"]]   featureCounts matrix
#     snakemake@input[["samples"]]  sample sheet (tsv)
#     snakemake@output[["..."]]     output paths
#     snakemake@params[["..."]]     settings from config.yaml
#     snakemake@log[[1]]            log file
#  We never touch commandArgs() — Snakemake wires it all up for us.
# =====================================================================

# ---- 1. Send all console + message output to the rule's log file ----
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(pheatmap)
    library(RColorBrewer)
})

# ---- 2. Pull parameters out of the snakemake object -----------------
condition_col   <- snakemake@params[["condition_col"]]
reference_level <- snakemake@params[["reference_level"]]
treatment_level <- snakemake@params[["treatment_level"]]
padj_cut        <- as.numeric(snakemake@params[["padj"]])
lfc_cut         <- as.numeric(snakemake@params[["lfc"]])
sample_order    <- snakemake@params[["sample_order"]]   # character vector

message("Contrast: ", treatment_level, " vs ", reference_level,
        "  (reference = ", reference_level, ")")
message("Thresholds: padj < ", padj_cut, ", |log2FC| >= ", lfc_cut)

# ---- 3. Read the featureCounts matrix -------------------------------
# featureCounts files start with a "# Program:..." comment line, so we
# skip it. Columns 1-6 are annotation (Geneid, Chr, Start, End, Strand,
# Length); the rest are one count column per BAM.
fc <- read.delim(snakemake@input[["counts"]], comment.char = "#",
                 check.names = FALSE)
rownames(fc) <- fc$Geneid
counts <- as.matrix(fc[, 7:ncol(fc), drop = FALSE])

# The columns are named by BAM path; rename them to sample names. The
# order matches because the rule passed the BAMs in `sample_order`.
stopifnot(ncol(counts) == length(sample_order))
colnames(counts) <- sample_order
mode(counts) <- "integer"

# ---- 4. Build colData (sample metadata) from the sample sheet -------
coldata <- read.delim(snakemake@input[["samples"]], check.names = FALSE)
rownames(coldata) <- coldata$sample
# Standardise the grouping column to "condition" so the design is simple.
coldata$condition <- as.character(coldata[[condition_col]])

# Keep only the two groups being compared (the sheet may hold more).
keep <- coldata$condition %in% c(reference_level, treatment_level)
coldata <- coldata[keep, , drop = FALSE]
# Make the reference level FIRST so log2FC = treatment vs reference.
coldata$condition <- factor(coldata$condition,
                            levels = c(reference_level, treatment_level))

# Align the count columns to the (subset of) samples, same order.
counts <- counts[, rownames(coldata), drop = FALSE]

message("Samples per group:")
print(table(coldata$condition))

# ---- 5. DESeqDataSet -------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ condition)
# Drop genes with zero counts everywhere — they carry no information.
dds <- dds[rowSums(counts(dds)) > 0, ]

# Do we have replicates in BOTH groups? DESeq2's statistics require >=2.
group_sizes <- table(coldata$condition)
has_replicates <- all(group_sizes >= 2)

# Helper: write a one-message PDF when a real plot can't be made.
placeholder_pdf <- function(path, msg) {
    pdf(path); plot.new()
    text(0.5, 0.5, msg, cex = 1.1); dev.off()
}

if (has_replicates) {
    # ============ FULL STATISTICAL ANALYSIS ==========================
    message("Replicates detected -> running full DESeq2 Wald test.")
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", treatment_level, reference_level),
                   alpha = padj_cut)
    res <- res[order(res$padj), ]

    res_df <- as.data.frame(res)
    res_df <- cbind(gene = rownames(res_df), res_df)

    # Significant = passes BOTH padj and |log2FC| cutoffs.
    sig <- subset(res_df, !is.na(padj) & padj < padj_cut &
                          abs(log2FoldChange) >= lfc_cut)

    # MA plot (mean expression vs log2FC).
    pdf(snakemake@output[["ma"]]); plotMA(res, ylim = c(-5, 5)); dev.off()

    # Volcano plot.
    vol <- res_df
    vol$signif <- with(vol, !is.na(padj) & padj < padj_cut &
                              abs(log2FoldChange) >= lfc_cut)
    p <- ggplot(vol, aes(log2FoldChange, -log10(padj))) +
        geom_point(aes(colour = signif), alpha = 0.6, size = 1) +
        scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                            name = "significant") +
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
        geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
        labs(title = paste(treatment_level, "vs", reference_level),
             x = "log2 fold change", y = "-log10 adjusted p-value") +
        theme_bw()
    ggsave(snakemake@output[["volcano"]], p, width = 7, height = 6)

} else {
    # ============ NO-REPLICATE FALLBACK ==============================
    # With 1 sample per group DESeq2 cannot estimate dispersion, so a
    # statistical test is impossible. We still give DESCRIPTIVE fold
    # changes from normalised counts so the pipeline completes — but
    # there are NO p-values, so treat these as exploratory only.
    message("WARNING: <2 replicates in a group. No statistics possible.")
    message("Producing descriptive log2 fold-changes ONLY (no p-values).")
    dds <- estimateSizeFactors(dds)
    norm <- counts(dds, normalized = TRUE)

    ref_samples   <- rownames(coldata)[coldata$condition == reference_level]
    treat_samples <- rownames(coldata)[coldata$condition == treatment_level]
    ref_mean   <- rowMeans(norm[, ref_samples, drop = FALSE])
    treat_mean <- rowMeans(norm[, treat_samples, drop = FALSE])

    res_df <- data.frame(
        gene           = rownames(norm),
        baseMean       = rowMeans(norm),
        log2FoldChange = log2((treat_mean + 1) / (ref_mean + 1)),
        lfcSE          = NA_real_,
        stat           = NA_real_,
        pvalue         = NA_real_,
        padj           = NA_real_,
        row.names      = NULL
    )
    res_df <- res_df[order(-abs(res_df$log2FoldChange)), ]

    # Without p-values we can only rank by effect size.
    sig <- subset(res_df, abs(log2FoldChange) >= lfc_cut)

    # MA-style plot: mean expression vs log2FC (no significance colour).
    p_ma <- ggplot(res_df, aes(log2(baseMean + 1), log2FoldChange)) +
        geom_point(alpha = 0.4, size = 1, colour = "grey40") +
        geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed") +
        labs(title = "MA plot (descriptive - no replicates)",
             x = "log2 mean normalised count", y = "log2 fold change") +
        theme_bw()
    ggsave(snakemake@output[["ma"]], p_ma, width = 7, height = 6)

    placeholder_pdf(snakemake@output[["volcano"]],
                    "Volcano plot needs p-values.\nNo replicates -> no statistics.")
}

# ---- 6. Write the tables (shared by both branches) ------------------
write.table(res_df, snakemake@output[["results"]],
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sig, snakemake@output[["significant"]],
            sep = "\t", quote = FALSE, row.names = FALSE)

norm_out <- counts(dds, normalized = TRUE)
norm_out <- cbind(gene = rownames(norm_out), as.data.frame(norm_out))
write.table(norm_out, snakemake@output[["normalized"]],
            sep = "\t", quote = FALSE, row.names = FALSE)

message("Significant genes: ", nrow(sig), " / ", nrow(res_df), " tested")

# ---- 7. Diagnostic plots (PCA + sample-distance heatmap) ------------
# We need a transform that puts counts on a comparable scale for PCA /
# distances. The variance-stabilising transform (vst) is best, but it
# estimates dispersion, which is impossible without replicates — there we
# fall back to a plain log2(normalised + 1) via normTransform().
if (has_replicates) {
    vsd <- tryCatch(
        vst(dds, blind = TRUE),
        error = function(e) varianceStabilizingTransformation(dds, blind = TRUE)
    )
} else {
    vsd <- normTransform(dds)   # log2(normalised counts + 1); no dispersion
}

# PCA - do samples separate by condition?
# geom_text (not geom_text_repel): ggrepel fails with "Viewport has zero
# dimension(s)" when there are only a couple of points to place, which is
# exactly the small-sample case this pipeline must survive.
p_pca <- plotPCA(vsd, intgroup = "condition") +
    geom_text(aes(label = name), size = 3, vjust = -0.6) +
    labs(title = "PCA of samples") + theme_bw()
ggsave(snakemake@output[["pca"]], p_pca, width = 7, height = 6)

# Sample-to-sample distance heatmap.
sampleDists <- dist(t(assay(vsd)))
distMat <- as.matrix(sampleDists)
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(distMat,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colours,
         main = "Sample-to-sample distances",
         filename = snakemake@output[["heatmap"]])

message("DESeq2 step complete.")
