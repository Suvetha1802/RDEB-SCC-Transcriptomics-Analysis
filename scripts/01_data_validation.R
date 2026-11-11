#!/usr/bin/env Rscript

# 01_data_validation.R
# Validate a gene-level raw count matrix and its sample metadata before DESeq2.
#
# Expected input files:
#   data/raw/gene_counts.tsv
#   data/metadata/sample_metadata.csv
#
# Count matrix format:
#   First column: gene_id
#   Remaining columns: one raw-count column per sample
#
# Metadata format (CSV):
#   sample_id,group
#   sample_01,Healthy
#   sample_02,RDEB
#   sample_03,RDEB_SCC
#
# Run from the repository root:
#   Rscript scripts/01_data_validation.R
# Or provide paths:
#   Rscript scripts/01_data_validation.R counts.tsv metadata.csv

options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
count_file <- if (length(args) >= 1) args[[1]] else "data/raw/gene_counts.tsv"
metadata_file <- if (length(args) >= 2) args[[2]] else "data/metadata/sample_metadata.csv"

expected_groups <- c("Healthy", "RDEB", "RDEB_SCC")
expected_group_sizes <- c(Healthy = 16L, RDEB = 10L, RDEB_SCC = 8L)
expected_total_samples <- sum(expected_group_sizes)

fail <- function(message) {
  stop(paste0("VALIDATION FAILED: ", message), call. = FALSE)
}

if (!file.exists(count_file)) {
  fail(paste0("Count matrix not found: ", count_file,
              ". Place the original raw gene-count matrix at this path or pass its path as argument 1."))
}
if (!file.exists(metadata_file)) {
  fail(paste0("Sample metadata not found: ", metadata_file,
              ". Create a CSV with columns sample_id and group, using the original sample IDs and group assignments."))
}

message("Reading count matrix: ", count_file)
counts <- tryCatch(
  read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE,
             quote = "", comment.char = ""),
  error = function(e) fail(paste("Could not read count matrix:", conditionMessage(e)))
)

message("Reading sample metadata: ", metadata_file)
metadata <- tryCatch(
  read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) fail(paste("Could not read metadata:", conditionMessage(e)))
)

if (ncol(counts) < 2L) {
  fail("Count matrix must contain a gene ID column and at least one sample column.")
}
if (!all(c("sample_id", "group") %in% names(metadata))) {
  fail("Metadata must contain columns named exactly 'sample_id' and 'group'.")
}
if (anyNA(metadata$sample_id) || any(trimws(metadata$sample_id) == "")) {
  fail("Metadata contains missing or blank sample_id values.")
}
if (anyDuplicated(metadata$sample_id)) {
  dup <- unique(metadata$sample_id[duplicated(metadata$sample_id)])
  fail(paste("Duplicate sample IDs in metadata:", paste(dup, collapse = ", ")))
}
if (anyNA(metadata$group) || any(trimws(metadata$group) == "")) {
  fail("Metadata contains missing or blank group values.")
}

metadata$sample_id <- trimws(as.character(metadata$sample_id))
metadata$group <- trimws(as.character(metadata$group))

unknown_groups <- setdiff(unique(metadata$group), expected_groups)
if (length(unknown_groups) > 0L) {
  fail(paste("Unexpected group label(s):", paste(unknown_groups, collapse = ", "),
             ". Use Healthy, RDEB, and RDEB_SCC consistently."))
}

group_counts <- table(factor(metadata$group, levels = expected_groups))
if (!identical(as.integer(group_counts), as.integer(expected_group_sizes))) {
  observed <- paste(names(group_counts), as.integer(group_counts), sep = "=", collapse = "; ")
  expected <- paste(names(expected_group_sizes), as.integer(expected_group_sizes),
                    sep = "=", collapse = "; ")
  fail(paste0("Group sample counts do not match the recorded study design. Observed: ",
              observed, ". Expected: ", expected, ". Confirm the original sample metadata; do not guess or duplicate samples."))
}

gene_id_col <- names(counts)[1]
gene_ids <- as.character(counts[[1]])

if (anyNA(gene_ids) || any(trimws(gene_ids) == "")) {
  fail("The first count-matrix column (gene ID) contains missing or blank values.")
}
if (anyDuplicated(gene_ids)) {
  fail("Duplicate gene IDs found in the first count-matrix column. Resolve duplicates using the original quantification output before DESeq2.")
}

count_sample_ids <- names(counts)[-1]
if (anyDuplicated(count_sample_ids)) {
  fail("Duplicate sample column names found in the count matrix.")
}
if (anyNA(count_sample_ids) || any(trimws(count_sample_ids) == "")) {
  fail("Count matrix contains missing or blank sample column names.")
}

missing_in_counts <- setdiff(metadata$sample_id, count_sample_ids)
missing_in_metadata <- setdiff(count_sample_ids, metadata$sample_id)

if (length(missing_in_counts) > 0L || length(missing_in_metadata) > 0L) {
  fail(paste0(
    "Sample IDs do not match between files. Missing from count matrix: ",
    ifelse(length(missing_in_counts), paste(missing_in_counts, collapse = ", "), "none"),
    ". Missing from metadata: ",
    ifelse(length(missing_in_metadata), paste(missing_in_metadata, collapse = ", "), "none"),
    ". Check exact sample names; do not rename based on assumptions."
  ))
}

if (length(count_sample_ids) != expected_total_samples) {
  fail(paste0("Expected ", expected_total_samples, " sample columns, but found ",
              length(count_sample_ids), "."))
}

count_values <- counts[-1]
numeric_columns <- vapply(count_values, is.numeric, logical(1))
if (!all(numeric_columns)) {
  fail(paste("Non-numeric count column(s):",
             paste(names(count_values)[!numeric_columns], collapse = ", ")))
}

count_matrix <- as.matrix(count_values)
if (anyNA(count_matrix) || any(!is.finite(count_matrix))) {
  fail("Count matrix contains NA, NaN, or infinite values.")
}
if (any(count_matrix < 0)) {
  fail("Count matrix contains negative values.")
}
if (any(count_matrix != floor(count_matrix))) {
  fail("Count matrix contains non-integer values. DESeq2 requires raw integer counts; do not round normalized expression values.")
}
if (any(colSums(count_matrix) == 0)) {
  fail("At least one sample has a total count of zero.")
}

# Reorder metadata to exactly match count-matrix sample-column order.
metadata <- metadata[match(count_sample_ids, metadata$sample_id), , drop = FALSE]
rownames(metadata) <- metadata$sample_id

report_dir <- "results/QC"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
report_file <- file.path(report_dir, "data_validation_report.txt")

report <- c(
  "RDEB-SCC data validation report",
  paste("Date:", format(Sys.Date(), "%Y-%m-%d")),
  paste("Count matrix:", normalizePath(count_file, winslash = "/", mustWork = FALSE)),
  paste("Metadata:", normalizePath(metadata_file, winslash = "/", mustWork = FALSE)),
  paste("Gene ID column:", gene_id_col),
  paste("Number of genes:", nrow(counts)),
  paste("Number of samples:", length(count_sample_ids)),
  "",
  "Sample counts by group:",
  capture.output(print(group_counts)),
  "",
  "Sample order used (matches count-matrix columns):",
  paste(metadata$sample_id, metadata$group, sep = " : "),
  "",
  "Checks passed:",
  "- Unique, non-empty gene IDs",
  "- Unique sample IDs",
  "- Exact sample ID match between count matrix and metadata",
  "- Expected group labels and sample counts",
  "- Numeric, finite, non-negative integer raw counts",
  "- No zero-total samples",
  "",
  "Status: PASS. This validates input structure only; it does not run differential expression."
)

writeLines(report, report_file)
cat(paste(report, collapse = "\n"), "\n")
message("Validation report written to: ", report_file)
