# Summarize core neutrophil marker genes in existing mRNA DE results
#
# This script does not rerun differential expression. It extracts the
# independently defined core mouse neutrophil marker panel from the existing
# treatment + sex + age mRNA limma results for the three planned contrasts.
#
# Inputs:
#   results/neutrophil_analysis/gene_sets/core_mouse_neutrophil_marker_coverage.csv
#     Source-traceable marker panel and measured-gene coverage.
#   results/mrna/differential_expression/treatment_sex_age_limma/
#     transcript_cluster_MCAO1hr_vs_Sham_mrna.csv
#     transcript_cluster_MCAO3hr_vs_Sham_mrna.csv
#     transcript_cluster_MCAO3hr_vs_MCAO1hr_mrna.csv
#     Existing mRNA results from expression ~ treatment + sex + age.
#
# Outputs:
#   results/neutrophil_analysis/tables/core_neutrophil_marker_mrna_de_results.csv
#     One row per measured marker and contrast with DE statistics.
#   results/neutrophil_analysis/tables/core_neutrophil_marker_mrna_de_summary.txt
#     Marker counts, threshold summaries, and interpretation safeguards.

options(stringsAsFactors = FALSE)

project_root <- normalizePath(".", mustWork = TRUE)
panel_path <- file.path(
  project_root,
  "results/neutrophil_analysis/gene_sets/core_mouse_neutrophil_marker_coverage.csv"
)
de_dir <- file.path(
  project_root, "results/mrna/differential_expression/treatment_sex_age_limma"
)
output_dir <- file.path(project_root, "results/neutrophil_analysis/tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_file <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  TRUE
}
invisible(required_file(panel_path))

contrast_files <- c(
  MCAO1hr_vs_Sham = "transcript_cluster_MCAO1hr_vs_Sham_mrna.csv",
  MCAO3hr_vs_Sham = "transcript_cluster_MCAO3hr_vs_Sham_mrna.csv",
  MCAO3hr_vs_MCAO1hr = "transcript_cluster_MCAO3hr_vs_MCAO1hr_mrna.csv"
)
invisible(vapply(file.path(de_dir, contrast_files), required_file, logical(1)))

panel <- read.csv(panel_path, check.names = FALSE)
required_panel_columns <- c(
  "gene_symbol", "panel_tier", "evidence", "measured_in_mrna_matrix"
)
missing_panel_columns <- setdiff(required_panel_columns, names(panel))
if (length(missing_panel_columns) > 0) {
  stop("Marker coverage file is missing columns: ",
       paste(missing_panel_columns, collapse = ", "), call. = FALSE)
}
panel <- panel[panel$measured_in_mrna_matrix, ]
panel$gene_symbol <- trimws(as.character(panel$gene_symbol))
panel <- panel[!is.na(panel$gene_symbol) & panel$gene_symbol != "", ]
panel <- panel[!duplicated(panel$gene_symbol), ]

read_de <- function(path, contrast_name) {
  de <- read.csv(path, check.names = FALSE)
  required_de_columns <- c("SYMBOL", "logFC", "P.Value", "adj.P.Val")
  missing_de_columns <- setdiff(required_de_columns, names(de))
  if (length(missing_de_columns) > 0) {
    stop("DE file is missing columns: ", paste(missing_de_columns, collapse = ", "),
         " (", path, ")", call. = FALSE)
  }
  de$gene_symbol <- trimws(as.character(de$SYMBOL))
  de <- de[!is.na(de$gene_symbol) & de$gene_symbol != "", ]
  de <- de[!duplicated(de$gene_symbol), ]
  out <- merge(
    panel,
    de[, c("gene_symbol", "logFC", "P.Value", "adj.P.Val")],
    by = "gene_symbol", all.x = TRUE, sort = FALSE
  )
  out$contrast <- contrast_name
  out$nominal_p_lt_0.05 <- !is.na(out$P.Value) & out$P.Value < 0.05
  out$FDR_lt_0.10 <- !is.na(out$adj.P.Val) & out$adj.P.Val < 0.10
  out$FDR_lt_0.05 <- !is.na(out$adj.P.Val) & out$adj.P.Val < 0.05
  out[, c("contrast", "gene_symbol", "panel_tier", "evidence", "logFC",
          "P.Value", "adj.P.Val", "nominal_p_lt_0.05", "FDR_lt_0.10",
          "FDR_lt_0.05")]
}

results <- do.call(
  rbind,
  Map(read_de, file.path(de_dir, contrast_files), names(contrast_files))
)
results <- results[order(results$contrast, results$adj.P.Val, results$P.Value), ]
rownames(results) <- NULL

write.csv(
  results,
  file.path(output_dir, "core_neutrophil_marker_mrna_de_results.csv"),
  row.names = FALSE,
  quote = TRUE
)

summary_lines <- c(
  "Core neutrophil marker mRNA DE summary",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Model: expression ~ treatment + sex + age",
  paste("Measured markers evaluated:", length(unique(panel$gene_symbol))),
  "",
  capture.output({
    for (contrast_name in names(contrast_files)) {
      x <- results[results$contrast == contrast_name, ]
      cat(contrast_name, "\n", sep = "")
      cat("  nominal p < 0.05: ", sum(x$nominal_p_lt_0.05), "\n", sep = "")
      cat("  FDR < 0.10: ", sum(x$FDR_lt_0.10), "\n", sep = "")
      cat("  FDR < 0.05: ", sum(x$FDR_lt_0.05), "\n", sep = "")
    }
  }),
  "",
  "Interpretation: the marker panel was defined independently of these DE",
  "results. Nominal p-values are exploratory; FDR values account for multiple",
  "testing in the complete mRNA analysis. Bulk expression cannot distinguish",
  "neutrophil abundance from altered transcription within neutrophils."
)
writeLines(
  summary_lines,
  file.path(output_dir, "core_neutrophil_marker_mrna_de_summary.txt")
)

message("Core neutrophil marker mRNA DE summary complete.")
message("Rows written: ", nrow(results))
