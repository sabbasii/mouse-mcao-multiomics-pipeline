# Integrate core neutrophil-marker mRNA and miRNA DE evidence
#
# This script joins results that were already generated separately. It does
# not rerun DE, target prediction, correlation, or mediation analysis.
# Inverse direction is treated as supportive compatibility only.
#
# Inputs:
#   results/neutrophil_analysis/tables/core_neutrophil_marker_mrna_de_results.csv
#   results/neutrophil_analysis/tables/neutrophil_marker_mirna_de_results.csv
#
# Outputs:
#   results/neutrophil_analysis/tables/neutrophil_marker_integrated_pairs.csv
#     One row per measured miRNA–marker–contrast with both-layer DE statistics.
#   results/neutrophil_analysis/tables/neutrophil_marker_integrated_summary.txt
#     Counts of inverse and DE-supported pairs by contrast.

options(stringsAsFactors = FALSE)
project_root <- normalizePath(".", mustWork = TRUE)
mrna_path <- file.path(project_root, "results/neutrophil_analysis/tables/core_neutrophil_marker_mrna_de_results.csv")
mirna_path <- file.path(project_root, "results/neutrophil_analysis/tables/neutrophil_marker_mirna_de_results.csv")
output_dir <- file.path(project_root, "results/neutrophil_analysis/tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_file <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  TRUE
}
invisible(required_file(mrna_path))
invisible(required_file(mirna_path))

mrna <- read.csv(mrna_path, check.names = FALSE)
mirna <- read.csv(mirna_path, check.names = FALSE)
mrna_needed <- c("contrast", "gene_symbol", "logFC", "P.Value", "adj.P.Val", "panel_tier")
mirna_needed <- c("contrast", "mirna_probeset_id", "mirna_id", "target_gene_symbol", "evidence_sources", "logFC", "P.Value", "adj.P.Val")
if (length(setdiff(mrna_needed, names(mrna))) > 0) stop("mRNA results are missing required columns.", call. = FALSE)
if (length(setdiff(mirna_needed, names(mirna))) > 0) stop("miRNA results are missing required columns.", call. = FALSE)

mrna <- mrna[, c("contrast", "gene_symbol", "panel_tier", "logFC", "P.Value", "adj.P.Val")]
names(mrna)[names(mrna) %in% c("logFC", "P.Value", "adj.P.Val")] <- c("mrna_logFC", "mrna_p_value", "mrna_fdr")
mirna <- mirna[, c("contrast", "mirna_probeset_id", "mirna_id", "target_gene_symbol", "panel_tier", "evidence_sources", "logFC", "P.Value", "adj.P.Val")]
names(mirna)[names(mirna) %in% c("logFC", "P.Value", "adj.P.Val")] <- c("mirna_logFC", "mirna_p_value", "mirna_fdr")

integrated <- merge(
  mirna,
  mrna,
  by.x = c("contrast", "target_gene_symbol"),
  by.y = c("contrast", "gene_symbol"),
  all = FALSE,
  sort = FALSE
)
integrated$mirna_nominal_p_lt_0.05 <- !is.na(integrated$mirna_p_value) & integrated$mirna_p_value < 0.05
integrated$mrna_nominal_p_lt_0.05 <- !is.na(integrated$mrna_p_value) & integrated$mrna_p_value < 0.05
integrated$mirna_fdr_lt_0.10 <- !is.na(integrated$mirna_fdr) & integrated$mirna_fdr < 0.10
integrated$mrna_fdr_lt_0.10 <- !is.na(integrated$mrna_fdr) & integrated$mrna_fdr < 0.10
integrated$inverse_direction <- !is.na(integrated$mirna_logFC) & !is.na(integrated$mrna_logFC) &
  sign(integrated$mirna_logFC) != sign(integrated$mrna_logFC)
integrated$both_nominal_p_lt_0.05 <- integrated$mirna_nominal_p_lt_0.05 & integrated$mrna_nominal_p_lt_0.05
integrated$inverse_and_both_nominal <- integrated$inverse_direction & integrated$both_nominal_p_lt_0.05
integrated$inverse_and_any_fdr_lt_0.10 <- integrated$inverse_direction &
  (integrated$mirna_fdr_lt_0.10 | integrated$mrna_fdr_lt_0.10)
integrated <- integrated[order(integrated$contrast, integrated$inverse_and_both_nominal, integrated$mirna_p_value), ]
rownames(integrated) <- NULL

write.csv(integrated, file.path(output_dir, "neutrophil_marker_integrated_pairs.csv"), row.names = FALSE, quote = TRUE)

summary_lines <- c(
  "Integrated core neutrophil-marker miRNA–mRNA summary", "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Integrated measured pairs:", nrow(integrated)),
  "",
  capture.output({
    for (contrast_name in unique(integrated$contrast)) {
      x <- integrated[integrated$contrast == contrast_name, ]
      cat(contrast_name, "\n", sep = "")
      cat("  inverse-direction pairs: ", sum(x$inverse_direction), "\n", sep = "")
      cat("  both layers nominal p < 0.05: ", sum(x$both_nominal_p_lt_0.05), "\n", sep = "")
      cat("  inverse and both layers nominal p < 0.05: ", sum(x$inverse_and_both_nominal), "\n", sep = "")
      cat("  inverse with at least one layer FDR < 0.10: ", sum(x$inverse_and_any_fdr_lt_0.10), "\n", sep = "")
    }
  }),
  "",
  "Inverse direction is compatible with repressive miRNA regulation but is",
  "not evidence of direct regulation, causality, or mediation. Pair-level",
  "results remain exploratory and retain the original target-evidence source."
)
writeLines(summary_lines, file.path(output_dir, "neutrophil_marker_integrated_summary.txt"))
message("Integrated neutrophil marker-pair summary complete.")
message("Pairs written: ", nrow(integrated))
