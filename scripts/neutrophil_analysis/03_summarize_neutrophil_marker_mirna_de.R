# Summarize miRNAs targeting the validated neutrophil marker panel
#
# This script does not rerun differential expression. It keeps target evidence
# only for the measured core mouse neutrophil marker panel and joins those
# measured miRNAs to the existing miRNA DE results.
#
# Inputs:
#   results/neutrophil_analysis/gene_sets/core_mouse_neutrophil_marker_coverage.csv
#   results/multiomics/miRNA_target_evidence/analysis_ready/mirna_mrna_target_pairs.csv
#   results/mirna/differential_expression/rma_normalized_mirna/treatment_sex_age_limma/
#     dabg20_<contrast>_mirna.csv
#
# Outputs:
#   results/neutrophil_analysis/tables/neutrophil_marker_mirna_de_results.csv
#   results/neutrophil_analysis/tables/neutrophil_marker_mirna_de_summary.txt

options(stringsAsFactors = FALSE)
project_root <- normalizePath(".", mustWork = TRUE)
panel_path <- file.path(project_root, "results/neutrophil_analysis/gene_sets/core_mouse_neutrophil_marker_coverage.csv")
target_path <- file.path(project_root, "results/multiomics/miRNA_target_evidence/analysis_ready/mirna_mrna_target_pairs.csv")
de_dir <- file.path(project_root, "results/mirna/differential_expression/rma_normalized_mirna/treatment_sex_age_limma")
output_dir <- file.path(project_root, "results/neutrophil_analysis/tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_file <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  TRUE
}
invisible(required_file(panel_path))
invisible(required_file(target_path))
contrast_files <- c(
  MCAO1hr_vs_Sham = "dabg20_MCAO1hr_vs_Sham_mirna.csv",
  MCAO3hr_vs_Sham = "dabg20_MCAO3hr_vs_Sham_mirna.csv",
  MCAO3hr_vs_MCAO1hr = "dabg20_MCAO3hr_vs_MCAO1hr_mirna.csv"
)
invisible(vapply(file.path(de_dir, contrast_files), required_file, logical(1)))

panel <- read.csv(panel_path, check.names = FALSE)
targets <- read.csv(target_path, check.names = FALSE)
panel_needed <- c("gene_symbol", "panel_tier", "measured_in_mrna_matrix")
target_needed <- c("mirna_id", "target_gene_symbol", "evidence_sources", "mirna_probeset_id")
if (length(setdiff(panel_needed, names(panel))) > 0) stop("Panel file is missing required columns.", call. = FALSE)
if (length(setdiff(target_needed, names(targets))) > 0) stop("Target file is missing required columns.", call. = FALSE)

panel <- panel[panel$measured_in_mrna_matrix, c("gene_symbol", "panel_tier")]
panel$gene_symbol <- trimws(as.character(panel$gene_symbol))
panel <- panel[!duplicated(panel$gene_symbol), ]
targets$target_gene_symbol <- trimws(as.character(targets$target_gene_symbol))
targets$mirna_id <- trimws(as.character(targets$mirna_id))
targets$mirna_probeset_id <- trimws(as.character(targets$mirna_probeset_id))
targets <- targets[!is.na(targets$mirna_probeset_id) & targets$mirna_probeset_id != "", ]
targets <- merge(targets, panel, by.x = "target_gene_symbol", by.y = "gene_symbol", all = FALSE, sort = FALSE)
targets <- targets[!is.na(targets$mirna_id) & targets$mirna_id != "", ]

# Collapse duplicate database rows while preserving evidence provenance.
targets$pair_key <- paste(targets$mirna_probeset_id, targets$target_gene_symbol, sep = "__")
targets <- aggregate(
  cbind(mirna_id, evidence_sources, panel_tier) ~ pair_key + mirna_probeset_id + target_gene_symbol,
  data = targets,
  FUN = function(x) paste(sort(unique(as.character(x))), collapse = ";")
)

read_de <- function(path, contrast_name) {
  de <- read.csv(path, check.names = FALSE)
  needed <- c("ProbeSetName", "logFC", "P.Value", "adj.P.Val")
  missing <- setdiff(needed, names(de))
  if (length(missing) > 0) stop("miRNA DE file is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  de$mirna_probeset_id <- trimws(as.character(de$ProbeSetName))
  de <- de[!duplicated(de$mirna_probeset_id), ]
  out <- merge(targets, de[, c("mirna_probeset_id", "logFC", "P.Value", "adj.P.Val")], by = "mirna_probeset_id", all.x = TRUE, sort = FALSE)
  out$contrast <- contrast_name
  out$nominal_p_lt_0.05 <- !is.na(out$P.Value) & out$P.Value < 0.05
  out$FDR_lt_0.10 <- !is.na(out$adj.P.Val) & out$adj.P.Val < 0.10
  out$FDR_lt_0.05 <- !is.na(out$adj.P.Val) & out$adj.P.Val < 0.05
  out[, c("contrast", "mirna_probeset_id", "mirna_id", "target_gene_symbol", "panel_tier", "evidence_sources", "logFC", "P.Value", "adj.P.Val", "nominal_p_lt_0.05", "FDR_lt_0.10", "FDR_lt_0.05")]
}

results <- do.call(rbind, Map(read_de, file.path(de_dir, contrast_files), names(contrast_files)))
results <- results[!is.na(results$logFC), ]
results <- results[order(results$contrast, results$adj.P.Val, results$P.Value), ]
rownames(results) <- NULL
write.csv(results, file.path(output_dir, "neutrophil_marker_mirna_de_results.csv"), row.names = FALSE, quote = TRUE)

summary_lines <- c(
  "Core neutrophil marker-targeting miRNA DE summary", "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "miRNA model: expression ~ treatment + sex + age",
  paste("Measured marker-target rows with DE results:", nrow(results)),
  paste("Unique measured miRNAs:", length(unique(results$mirna_probeset_id))),
  paste("Unique marker genes targeted:", length(unique(results$target_gene_symbol))), "",
  capture.output({
    for (contrast_name in names(contrast_files)) {
      x <- results[results$contrast == contrast_name, ]
      cat(contrast_name, "\n", sep = "")
      cat("  pairs with nominal p < 0.05: ", sum(x$nominal_p_lt_0.05), "\n", sep = "")
      cat("  pairs with FDR < 0.10: ", sum(x$FDR_lt_0.10), "\n", sep = "")
      cat("  pairs with FDR < 0.05: ", sum(x$FDR_lt_0.05), "\n", sep = "")
      cat("  unique miRNAs with nominal p < 0.05: ", length(unique(x$mirna_probeset_id[x$nominal_p_lt_0.05])), "\n", sep = "")
    }
  }),
  "",
  "Target evidence is retained as provenance; miRTarBase and TargetScanMouse",
  "are not treated as equivalent evidence. This analysis does not establish",
  "direct miRNA regulation or causality."
)
writeLines(summary_lines, file.path(output_dir, "neutrophil_marker_mirna_de_summary.txt"))
message("Core neutrophil marker-targeting miRNA DE summary complete.")
message("Rows written: ", nrow(results))
