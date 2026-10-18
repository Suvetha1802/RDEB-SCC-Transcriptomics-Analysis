# RDEB-SCC transcriptomics: three-group DESeq2 analysis
# Requires: DESeq2, ggplot2
# Input counts: rows = genes, columns = samples, first column = gene identifier.
# Input metadata: columns sample_id, condition; condition must be Healthy, RDEB, or SCC.
# IMPORTANT: replace sample IDs in metadata with exact count-matrix column names.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# ---------------- USER SETTINGS ----------------
count_file <- "D:/final rdeb - scc/counts/merged_counts_all_samples.csv"
metadata_file <- "D:/final rdeb - scc/counts/sample_metadata.csv"
outdir <- "D:/final rdeb - scc/DESeq2_corrected"
padj_cutoff <- 0.05
lfc_cutoff <- 1
min_total_count <- 10
# ------------------------------------------------

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(count_file)) stop("Count matrix not found: ", count_file)
if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

raw_counts <- read.csv(count_file, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(raw_counts) < 3) stop("Count matrix must contain a gene-ID column and at least two sample columns.")
gene_col <- names(raw_counts)[1]
gene_ids <- trimws(as.character(raw_counts[[gene_col]]))
if (anyNA(gene_ids) || any(gene_ids == "")) stop("Gene identifiers contain missing/blank values.")
if (anyDuplicated(gene_ids)) stop("Duplicate gene identifiers found. Resolve duplicates before DESeq2: ", paste(unique(gene_ids[duplicated(gene_ids)]), collapse = ", "))

count_df <- raw_counts[, -1, drop = FALSE]
if (any(!vapply(count_df, is.numeric, logical(1)))) {
  converted <- lapply(count_df, function(x) suppressWarnings(as.numeric(as.character(x))))
  bad <- names(converted)[vapply(converted, function(x) anyNA(x), logical(1))]
  if (length(bad)) stop("Non-numeric or missing counts in sample columns: ", paste(bad, collapse = ", "))
  count_df <- as.data.frame(converted, check.names = FALSE)
}
count_mat <- as.matrix(count_df)
rownames(count_mat) <- gene_ids
if (any(!is.finite(count_mat)) || any(count_mat < 0)) stop("Counts must be finite and non-negative.")
if (any(abs(count_mat - round(count_mat)) > 1e-6)) stop("DESeq2 requires raw integer counts, not TPM/FPKM/normalized counts.")
storage.mode(count_mat) <- "integer"

meta <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
required_cols <- c("sample_id", "condition")
if (!all(required_cols %in% names(meta))) stop("Metadata must have columns named sample_id and condition.")
meta$sample_id <- trimws(as.character(meta$sample_id))
meta$condition <- trimws(as.character(meta$condition))
if (anyDuplicated(meta$sample_id)) stop("Duplicate sample_id values in metadata.")
if (!setequal(colnames(count_mat), meta$sample_id)) {
  missing_meta <- setdiff(colnames(count_mat), meta$sample_id)
  extra_meta <- setdiff(meta$sample_id, colnames(count_mat))
  stop("Count columns and metadata sample IDs do not match. Missing metadata IDs: ",
       paste(missing_meta, collapse = ", "), "; metadata IDs absent from counts: ",
       paste(extra_meta, collapse = ", "))
}
meta <- meta[match(colnames(count_mat), meta$sample_id), , drop = FALSE]
stopifnot(identical(meta$sample_id, colnames(count_mat)))
expected_groups <- c("Healthy", "RDEB", "SCC")
if (!all(meta$condition %in% expected_groups)) stop("condition values must be exactly: Healthy, RDEB, SCC.")
meta$condition <- factor(meta$condition, levels = expected_groups)
observed_n <- table(meta$condition)
message("Sample counts by group:\n", paste(capture.output(print(observed_n)), collapse = "\n"))
expected_n <- c(Healthy = 16, RDEB = 10, SCC = 8)
if (!identical(as.integer(observed_n[names(expected_n)]), as.integer(expected_n))) {
  stop("Group sizes do not match expected design (Healthy=16, RDEB=10, SCC=8). Correct sample_metadata.csv if this is the intended dataset.")
}

write.csv(meta, file.path(outdir, "Validated_Sample_Metadata.csv"), row.names = FALSE)
keep <- rowSums(count_mat) >= min_total_count
count_mat <- count_mat[keep, , drop = FALSE]
if (nrow(count_mat) < 1) stop("No genes remain after count filtering.")

coldata <- as.data.frame(meta[, c("sample_id", "condition")])
rownames(coldata) <- coldata$sample_id
coldata$sample_id <- NULL
dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = coldata, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "Healthy")
dds <- DESeq(dds)
saveRDS(dds, file.path(outdir, "DESeq2_object.rds"))

# PCA from variance-stabilized counts
vsd <- vst(dds, blind = FALSE)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))
p_pca <- ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  labs(title = "RDEB-SCC PCA", x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) + theme_bw()
ggsave(file.path(outdir, "PCA_all_samples.png"), p_pca, width = 8, height = 6, dpi = 300)

# Direction means the first group listed is compared with the second group.
contrasts <- list(
  RDEB_vs_Healthy = c("RDEB", "Healthy"),
  SCC_vs_Healthy = c("SCC", "Healthy"),
  SCC_vs_RDEB = c("SCC", "RDEB")
)

for (contrast_name in names(contrasts)) {
  pair <- contrasts[[contrast_name]]
  res <- results(dds, contrast = c("condition", pair[1], pair[2]), alpha = padj_cutoff)
  res_df <- as.data.frame(res)
  res_df$GeneID <- rownames(res_df)
  res_df <- res_df[, c("GeneID", setdiff(names(res_df), "GeneID")), drop = FALSE]
  res_df$Significance <- "Not significant"
  res_df$Significance[!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >= lfc_cutoff] <- "Up"
  res_df$Significance[!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange <= -lfc_cutoff] <- "Down"

  sig <- res_df[!is.na(res_df$padj) & res_df$padj < padj_cutoff & abs(res_df$log2FoldChange) >= lfc_cutoff, , drop = FALSE]
  up <- sig[sig$log2FoldChange >= lfc_cutoff, , drop = FALSE]
  down <- sig[sig$log2FoldChange <= -lfc_cutoff, , drop = FALSE]
  cdir <- file.path(outdir, contrast_name)
  dir.create(cdir, showWarnings = FALSE)
  write.csv(res_df, file.path(cdir, "All_DESeq2_Results.csv"), row.names = FALSE)
  write.csv(sig, file.path(cdir, "Significant_DEGs.csv"), row.names = FALSE)
  write.csv(up, file.path(cdir, "Upregulated_DEGs.csv"), row.names = FALSE)
  write.csv(down, file.path(cdir, "Downregulated_DEGs.csv"), row.names = FALSE)
  write.csv(head(up[order(-up$log2FoldChange), , drop = FALSE], 50), file.path(cdir, "Top50_Upregulated.csv"), row.names = FALSE)
  write.csv(head(down[order(down$log2FoldChange), , drop = FALSE], 50), file.path(cdir, "Top50_Downregulated.csv"), row.names = FALSE)

  pdf(file.path(cdir, "MA_plot.pdf"), width = 7, height = 5)
  plotMA(res, alpha = padj_cutoff, main = contrast_name)
  dev.off()

  volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pmax(padj, 1e-300)), color = Significance)) +
    geom_point(alpha = 0.65, size = 1.2, na.rm = TRUE) +
    scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "Not significant" = "grey70")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = contrast_name, x = "log2 fold change", y = "-log10 adjusted p-value") + theme_bw()
  ggsave(file.path(cdir, "Volcano_plot.png"), volcano, width = 7, height = 6, dpi = 300)
  message(contrast_name, ": ", nrow(sig), " significant DEGs; ", nrow(up), " up; ", nrow(down), " down.")
}

# Save normalized counts for downstream visualization/integration (not as DESeq2 input).
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
norm_counts$GeneID <- rownames(norm_counts)
norm_counts <- norm_counts[, c("GeneID", setdiff(names(norm_counts), "GeneID")), drop = FALSE]
write.csv(norm_counts, file.path(outdir, "Normalized_Counts_for_Visualization.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Finished. Review QC/PCA and group assignments before interpreting DEGs.")
