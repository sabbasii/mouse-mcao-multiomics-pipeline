#!/usr/bin/env Rscript

# Summarize treatment-associated miRNA--mRNA target pairs. This script adds
# transparent FDR and direction flags; it does not perform biological
# interpretation or pathway analysis.
#
# Input:
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_de_supported_target_pairs.csv
#
# Output:
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_de_pair_evidence_summary.csv
#
# Usage:
#   Rscript scripts/multiomics/06_summarize_target_pair_de_results.R

input_file <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "analysis_ready",
  "mirna_mrna_de_supported_target_pairs.csv"
)
output_dir <- dirname(input_file)
output_file <- file.path(
  output_dir, "mirna_mrna_de_pair_evidence_summary.csv"
)

if (!file.exists(input_file)) {
  stop("Missing DE-supported pair table: ", input_file, call. = FALSE)
}

pairs <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)

required_columns <- c(
  "contrast", "mirna_logFC", "mirna_adj.P.Val",
  "mrna_logFC", "mrna_adj.P.Val"
)
missing_columns <- setdiff(required_columns, names(pairs))
if (length(missing_columns) > 0L) {
  stop("Input is missing column(s): ", paste(missing_columns, collapse = ", "),
       call. = FALSE)
}

fdr_cutoff <- 0.10
pairs$mirna_fdr_lt_0_10 <- !is.na(pairs$mirna_adj.P.Val) &
  pairs$mirna_adj.P.Val < fdr_cutoff
pairs$mrna_fdr_lt_0_10 <- !is.na(pairs$mrna_adj.P.Val) &
  pairs$mrna_adj.P.Val < fdr_cutoff
pairs$both_fdr_lt_0_10 <- pairs$mirna_fdr_lt_0_10 & pairs$mrna_fdr_lt_0_10

pairs$inverse_logFC_direction <- !is.na(pairs$mirna_logFC) &
  !is.na(pairs$mrna_logFC) &
  pairs$mirna_logFC != 0 & pairs$mrna_logFC != 0 &
  sign(pairs$mirna_logFC) != sign(pairs$mrna_logFC)

pairs$candidate_category <- ifelse(
  pairs$both_fdr_lt_0_10,
  "both_layers_fdr_lt_0.10",
  ifelse(pairs$mirna_fdr_lt_0_10, "mirna_only_fdr_lt_0.10",
    ifelse(pairs$mrna_fdr_lt_0_10, "mrna_only_fdr_lt_0.10",
      "neither_layer_fdr_lt_0.10"))
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(pairs, output_file, row.names = FALSE, na = "")

cat("FDR cutoff: ", fdr_cutoff, "\n", sep = "")
cat("Wrote ", nrow(pairs), " pair-summary rows: ", output_file, "\n", sep = "")
cat("\nCandidate categories by contrast:\n")
print(table(pairs$contrast, pairs$candidate_category))
cat("\nInverse-direction pairs by contrast:\n")
print(with(pairs, table(contrast, inverse_logFC_direction)))
