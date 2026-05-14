#!/usr/bin/env Rscript

# Run sample-level quality control on the primary analysis-ready mRNA matrix.
#
# This script reviews normalized expression distributions, PCA, unsupervised
# clustering, sample-to-sample correlation, and sample-level distance metrics.
# Statistical flags identify samples for manual review only. The script never
# excludes a sample automatically.
#
# Prerequisite:
#   Rscript scripts/mrna/03_prepare_analysis_ready_mrna.R
#
# Inputs:
#   results/mrna/analysis_ready/expression_matrix_unique_gene_mapped_mrna.csv
#   results/mrna/analysis_ready/analysis_samples_mrna.csv
#
# Outputs:
#   results/mrna/qc/analysis_ready/
#     expression_boxplot_mrna.png
#     expression_density_mrna.png
#     pca_by_treatment_mrna.png
#     pca_by_sex_mrna.png
#     pca_by_age_group_mrna.png
#     pca_by_array_prefix_mrna.png
#     sample_clustering_mrna.png
#     pca_scores_mrna.csv
#     sample_qc_metrics_mrna.csv
#     possible_outliers_mrna.csv
#     sample_qc_summary_mrna.txt
#
# Run from the repository root:
#   Rscript scripts/mrna/04_sample_qc_mrna.R

expression_file <- file.path(
  "results", "mrna", "analysis_ready",
  "expression_matrix_unique_gene_mapped_mrna.csv"
)
sample_file <- file.path(
  "results", "mrna", "analysis_ready",
  "analysis_samples_mrna.csv"
)
output_dir <- file.path(
  "results", "mrna", "qc", "analysis_ready"
)

for (path in c(expression_file, sample_file)) {
  if (!file.exists(path)) {
    stop("Required input not found: ", path, call. = FALSE)
  }
}

expression_table <- read.csv(
  expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
samples <- read.csv(
  sample_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)

if (!"transcript_cluster_id" %in% names(expression_table)) {
  stop("Expression matrix lacks transcript_cluster_id.", call. = FALSE)
}
required_sample_columns <- c(
  "file_name",
  "animal_id",
  "sex",
  "age_group",
  "treatment",
  "is_control",
  "match_status"
)
missing_sample_columns <- setdiff(
  required_sample_columns,
  names(samples)
)
if (length(missing_sample_columns) > 0L) {
  stop(
    "Sample metadata lacks required fields: ",
    paste(missing_sample_columns, collapse = ", "),
    call. = FALSE
  )
}
if (any(samples$is_control) || any(samples$match_status != "matched")) {
  stop(
    "Analysis sample table contains a control or unmatched sample.",
    call. = FALSE
  )
}
if (anyDuplicated(samples$file_name)) {
  stop("Analysis sample filenames are duplicated.", call. = FALSE)
}
if (!all(samples$file_name %in% names(expression_table))) {
  stop(
    "Not every analysis sample is present in the expression matrix.",
    call. = FALSE
  )
}

expression <- as.matrix(
  expression_table[, samples$file_name, drop = FALSE]
)
storage.mode(expression) <- "numeric"
rownames(expression) <- as.character(
  expression_table$transcript_cluster_id
)
if (any(!is.finite(expression))) {
  stop("Expression matrix contains non-finite values.", call. = FALSE)
}
if (anyDuplicated(rownames(expression))) {
  stop("Expression matrix contains duplicated transcript clusters.", call. = FALSE)
}

extract_array_prefix <- function(file_names) {
  prefix <- sub(
    "^.*_([0-9]+)_[A-H][0-9]+[.]CEL$",
    "\\1",
    file_names,
    ignore.case = TRUE
  )
  prefix[prefix == file_names] <- NA_character_
  prefix
}
samples$array_prefix <- extract_array_prefix(samples$file_name)

flag_three_sd <- function(values, direction = c("high", "low")) {
  direction <- match.arg(direction)
  center <- mean(values, na.rm = TRUE)
  spread <- stats::sd(values, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) {
    return(rep(FALSE, length(values)))
  }
  if (direction == "high") {
    values > center + 3 * spread
  } else {
    values < center - 3 * spread
  }
}

pca <- stats::prcomp(
  t(expression),
  center = TRUE,
  scale. = FALSE
)
pca_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_scores <- cbind(
  samples,
  data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC1_percent = pca_variance[[1]],
    PC2_percent = pca_variance[[2]],
    row.names = NULL,
    check.names = FALSE
  )
)
pca_scores$pca_distance <- sqrt(
  pca_scores$PC1^2 + pca_scores$PC2^2
)

sample_correlations <- stats::cor(
  expression,
  use = "pairwise.complete.obs"
)
diag(sample_correlations) <- NA_real_

sample_metrics <- cbind(
  samples,
  data.frame(
    expression_mean = colMeans(expression),
    expression_median = apply(expression, 2, stats::median),
    expression_sd = apply(expression, 2, stats::sd),
    expression_IQR = apply(expression, 2, stats::IQR),
    distance_from_centroid = sqrt(
      colSums((expression - rowMeans(expression))^2)
    ),
    median_sample_correlation = apply(
      sample_correlations,
      2,
      stats::median,
      na.rm = TRUE
    ),
    pca_distance = pca_scores$pca_distance,
    row.names = NULL,
    check.names = FALSE
  )
)
sample_metrics$flag_high_expression_distance <- flag_three_sd(
  sample_metrics$distance_from_centroid,
  "high"
)
sample_metrics$flag_high_pca_distance <- flag_three_sd(
  sample_metrics$pca_distance,
  "high"
)
sample_metrics$flag_low_sample_correlation <- flag_three_sd(
  sample_metrics$median_sample_correlation,
  "low"
)
flag_columns <- grep(
  "^flag_",
  names(sample_metrics),
  value = TRUE
)
sample_metrics$possible_outlier_review <- rowSums(
  sample_metrics[, flag_columns, drop = FALSE]
) > 0L

plot_pca <- function(variable, path, display_name) {
  group <- factor(pca_scores[[variable]])
  palette <- grDevices::hcl.colors(
    max(3L, nlevels(group)),
    "Dark 3"
  )
  colors <- palette[as.integer(group)]

  grDevices::png(
    path,
    width = 10,
    height = 7.5,
    units = "in",
    res = 200
  )
  graphics::plot(
    pca_scores$PC1,
    pca_scores$PC2,
    col = colors,
    pch = 19,
    xlab = paste0(
      "PC1 (",
      round(pca_variance[[1]], 2),
      "%)"
    ),
    ylab = paste0(
      "PC2 (",
      round(pca_variance[[2]], 2),
      "%)"
    ),
    main = paste("Analysis-ready mRNA PCA by", display_name)
  )
  graphics::text(
    pca_scores$PC1,
    pca_scores$PC2,
    labels = pca_scores$animal_id,
    pos = 3,
    cex = 0.55
  )
  graphics::legend(
    "topright",
    legend = levels(group),
    col = palette[seq_len(nlevels(group))],
    pch = 19,
    cex = 0.8
  )
  grDevices::dev.off()
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

grDevices::png(
  file.path(output_dir, "expression_boxplot_mrna.png"),
  width = 14,
  height = 8,
  units = "in",
  res = 200
)
graphics::boxplot(
  expression,
  names = samples$animal_id,
  las = 2,
  outline = FALSE,
  cex.axis = 0.65,
  main = "Analysis-ready mRNA expression distributions",
  ylab = "RMA expression (log2)"
)
grDevices::dev.off()

density_colors <- grDevices::adjustcolor(
  grDevices::hcl.colors(ncol(expression), "Dynamic"),
  alpha.f = 0.65
)
grDevices::png(
  file.path(output_dir, "expression_density_mrna.png"),
  width = 10,
  height = 7,
  units = "in",
  res = 200
)
first_density <- stats::density(expression[, 1])
graphics::plot(
  first_density,
  col = density_colors[[1]],
  lwd = 1,
  main = "Analysis-ready mRNA expression densities",
  xlab = "RMA expression (log2)"
)
for (index in 2:ncol(expression)) {
  graphics::lines(
    stats::density(expression[, index]),
    col = density_colors[[index]],
    lwd = 1
  )
}
grDevices::dev.off()

plot_pca(
  "treatment",
  file.path(output_dir, "pca_by_treatment_mrna.png"),
  "treatment"
)
plot_pca(
  "sex",
  file.path(output_dir, "pca_by_sex_mrna.png"),
  "sex"
)
plot_pca(
  "age_group",
  file.path(output_dir, "pca_by_age_group_mrna.png"),
  "age group"
)
plot_pca(
  "array_prefix",
  file.path(output_dir, "pca_by_array_prefix_mrna.png"),
  "array prefix"
)

grDevices::png(
  file.path(output_dir, "sample_clustering_mrna.png"),
  width = 11,
  height = 9,
  units = "in",
  res = 200
)
graphics::plot(
  stats::hclust(stats::dist(t(expression))),
  labels = samples$animal_id,
  main = "Analysis-ready mRNA sample clustering",
  xlab = "",
  sub = ""
)
grDevices::dev.off()

write.csv(
  pca_scores,
  file.path(output_dir, "pca_scores_mrna.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
write.csv(
  sample_metrics,
  file.path(output_dir, "sample_qc_metrics_mrna.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
write.csv(
  sample_metrics[
    sample_metrics$possible_outlier_review,
    ,
    drop = FALSE
  ],
  file.path(output_dir, "possible_outliers_mrna.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

flagged_ids <- sample_metrics$animal_id[
  sample_metrics$possible_outlier_review
]
writeLines(
  c(
    "Dataset: analysis-ready Clariom S Mouse HT mRNA",
    paste("Expression input:", expression_file),
    paste("Sample input:", sample_file),
    paste(
      "QC expression dimensions:",
      paste(dim(expression), collapse = " x ")
    ),
    paste("Missing values:", sum(is.na(expression))),
    paste("NaN values:", sum(is.nan(expression))),
    paste("Infinite values:", sum(is.infinite(expression))),
    paste("Transcript clusters used:", nrow(expression)),
    paste("Biological samples used:", ncol(expression)),
    paste("PC1 variance (%):", round(pca_variance[[1]], 2)),
    paste("PC2 variance (%):", round(pca_variance[[2]], 2)),
    paste("Samples flagged for review:", length(flagged_ids)),
    paste(
      "Flagged animal IDs:",
      if (length(flagged_ids) == 0L) {
        "none"
      } else {
        paste(flagged_ids, collapse = ", ")
      }
    ),
    "Flag rule: more than 3 SD for expression distance, PCA distance, or low median sample correlation.",
    "No samples are excluded automatically by this QC script."
  ),
  file.path(output_dir, "sample_qc_summary_mrna.txt")
)

message("Completed analysis-ready mRNA sample QC.")
message("QC output: ", output_dir)
