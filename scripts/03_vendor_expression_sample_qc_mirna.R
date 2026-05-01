#!/usr/bin/env Rscript

# Sample-level QC annotation for miRNA 150001 vendor RMA-DABG expression.
# This script flags possible outliers for review only. It does not exclude
# samples automatically and does not run differential expression.

expression_file <- file.path(
  "results",
  "normalized",
  "vendor_rma_expression_mirna_150001.csv"
)
dabg_file <- file.path(
  "results",
  "normalized",
  "vendor_dabg_pvalues_mirna_150001.csv"
)
sample_sheet_file <- file.path(
  "results",
  "sample_sheet",
  "sample_sheet_mirna_150001.csv"
)
qc_root_dir <- file.path("results", "qc")
qc_dir <- file.path(qc_root_dir, "vendor_sample_qc")

pca_treatment_file <- file.path(
  qc_dir,
  "vendor_rma_pca_by_treatment_mirna_150001.png"
)
pca_age_file <- file.path(qc_dir, "vendor_rma_pca_by_age_mirna_150001.png")
pca_sex_file <- file.path(qc_dir, "vendor_rma_pca_by_sex_mirna_150001.png")
pca_array_prefix_file <- file.path(
  qc_dir,
  "vendor_rma_pca_by_array_prefix_mirna_150001.png"
)
pca_scores_file <- file.path(
  qc_dir,
  "vendor_rma_pca_scores_annotated_mirna_150001.csv"
)
clustering_file <- file.path(
  qc_dir,
  "vendor_rma_sample_clustering_annotated_mirna_150001.png"
)
sample_qc_file <- file.path(
  qc_dir,
  "vendor_rma_sample_qc_metrics_mirna_150001.csv"
)
outlier_file <- file.path(
  qc_dir,
  "vendor_rma_possible_outliers_mirna_150001.csv"
)
summary_file <- file.path(
  qc_dir,
  "vendor_rma_sample_qc_summary_mirna_150001.txt"
)

required_sample_columns <- c(
  "file_name",
  "animal_id",
  "sex",
  "age",
  "treatment",
  "is_control",
  "match_status"
)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

read_probe_matrix <- function(path, label) {
  table <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  if (!"ProbeSetName" %in% names(table)) {
    stop(label, " is missing ProbeSetName column: ", path, call. = FALSE)
  }

  probe_ids <- table$ProbeSetName
  matrix_values <- as.matrix(table[, setdiff(names(table), "ProbeSetName")])
  storage.mode(matrix_values) <- "numeric"
  rownames(matrix_values) <- probe_ids
  matrix_values
}

extract_array_prefix <- function(file_name) {
  sub("_.*$", "", sub("[.]CEL$", "", file_name, ignore.case = TRUE))
}

extract_array_position <- function(file_name) {
  stem <- sub("[.]CEL$", "", file_name, ignore.case = TRUE)
  matches <- regmatches(stem, regexpr("_[A-H][0-9]{2}_", stem))
  ifelse(
    nzchar(matches),
    gsub("_", "", matches),
    NA_character_
  )
}

flag_three_sd <- function(values, high = TRUE, low = FALSE) {
  values <- as.numeric(values)
  center <- mean(values, na.rm = TRUE)
  spread <- stats::sd(values, na.rm = TRUE)

  if (!is.finite(spread) || spread == 0) {
    return(rep(FALSE, length(values)))
  }

  flags <- rep(FALSE, length(values))
  if (high) {
    flags <- flags | values > center + 3 * spread
  }
  if (low) {
    flags <- flags | values < center - 3 * spread
  }
  flags
}

plot_pca_by_variable <- function(pca_scores, color_column, output_file, title) {
  groups <- factor(pca_scores[[color_column]])
  colors <- seq_along(levels(groups))

  png(output_file, width = 8, height = 7, units = "in", res = 300)
  plot(
    pca_scores$PC1,
    pca_scores$PC2,
    pch = 19,
    col = colors[groups],
    xlab = paste0("PC1 (", round(unique(pca_scores$PC1_percent), 1), "%)"),
    ylab = paste0("PC2 (", round(unique(pca_scores$PC2_percent), 1), "%)"),
    main = title
  )
  legend(
    "topright",
    legend = levels(groups),
    col = colors,
    pch = 19,
    cex = 0.8
  )
  dev.off()
}

stop_if_missing(expression_file, "Vendor RMA expression matrix")
stop_if_missing(dabg_file, "Vendor DABG p-value matrix")
stop_if_missing(sample_sheet_file, "Sample sheet")

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

expression_matrix <- read_probe_matrix(expression_file, "Vendor RMA expression matrix")
dabg_matrix <- read_probe_matrix(dabg_file, "Vendor DABG p-value matrix")

if (!identical(rownames(expression_matrix), rownames(dabg_matrix))) {
  stop("Expression and DABG matrices have different ProbeSetName rows.", call. = FALSE)
}

if (!identical(colnames(expression_matrix), colnames(dabg_matrix))) {
  stop("Expression and DABG matrices have different sample columns.", call. = FALSE)
}

if (any(!is.finite(expression_matrix))) {
  stop("Expression matrix contains non-finite values.", call. = FALSE)
}

if (any(!is.finite(dabg_matrix))) {
  stop("DABG p-value matrix contains non-finite values.", call. = FALSE)
}

sample_sheet <- read.csv(
  sample_sheet_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_sample_columns <- setdiff(required_sample_columns, names(sample_sheet))
if (length(missing_sample_columns) > 0) {
  stop(
    "Sample sheet is missing required columns: ",
    paste(missing_sample_columns, collapse = ", "),
    call. = FALSE
  )
}

sample_metadata <- sample_sheet[
  sample_sheet$is_control == FALSE &
    sample_sheet$match_status == "matched",
]

missing_metadata <- setdiff(colnames(expression_matrix), sample_metadata$file_name)
if (length(missing_metadata) > 0) {
  stop(
    "Expression matrix samples are missing from sample sheet:\n",
    paste(missing_metadata, collapse = "\n"),
    call. = FALSE
  )
}

extra_metadata <- setdiff(sample_metadata$file_name, colnames(expression_matrix))
if (length(extra_metadata) > 0) {
  stop(
    "Sample sheet non-control samples are missing from expression matrix:\n",
    paste(extra_metadata, collapse = "\n"),
    call. = FALSE
  )
}

sample_metadata <- sample_metadata[
  match(colnames(expression_matrix), sample_metadata$file_name),
]
rownames(sample_metadata) <- sample_metadata$file_name
sample_metadata$array_prefix <- extract_array_prefix(sample_metadata$file_name)
sample_metadata$array_position <- extract_array_position(sample_metadata$file_name)

pca <- prcomp(t(expression_matrix), center = TRUE, scale. = FALSE)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)

pca_scores <- data.frame(
  file_name = rownames(pca$x),
  animal_id = sample_metadata$animal_id,
  treatment = sample_metadata$treatment,
  age = sample_metadata$age,
  sex = sample_metadata$sex,
  array_prefix = sample_metadata$array_prefix,
  array_position = sample_metadata$array_position,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  PC1_percent = 100 * pca_var[1],
  PC2_percent = 100 * pca_var[2],
  PC3_percent = 100 * pca_var[3],
  stringsAsFactors = FALSE
)

write.csv(pca_scores, pca_scores_file, row.names = FALSE, quote = TRUE)

plot_pca_by_variable(
  pca_scores,
  "treatment",
  pca_treatment_file,
  "miRNA 150001 vendor RMA PCA by treatment"
)
plot_pca_by_variable(
  pca_scores,
  "age",
  pca_age_file,
  "miRNA 150001 vendor RMA PCA by age"
)
plot_pca_by_variable(
  pca_scores,
  "sex",
  pca_sex_file,
  "miRNA 150001 vendor RMA PCA by sex"
)
plot_pca_by_variable(
  pca_scores,
  "array_prefix",
  pca_array_prefix_file,
  "miRNA 150001 vendor RMA PCA by array prefix"
)

sample_dist <- dist(t(expression_matrix))
sample_hclust <- hclust(sample_dist, method = "complete")
cluster_labels <- paste(
  sample_metadata$animal_id,
  sample_metadata$treatment,
  sep = " | "
)

png(clustering_file, width = 12, height = 8, units = "in", res = 300)
plot(
  sample_hclust,
  labels = cluster_labels,
  main = "miRNA 150001 vendor RMA sample clustering",
  xlab = "",
  sub = "",
  cex = 0.7
)
dev.off()

expression_summaries <- t(apply(
  expression_matrix,
  2,
  function(values) {
    c(
      expr_min = min(values),
      expr_q1 = unname(stats::quantile(values, 0.25)),
      expr_median = median(values),
      expr_mean = mean(values),
      expr_q3 = unname(stats::quantile(values, 0.75)),
      expr_max = max(values),
      expr_sd = stats::sd(values),
      expr_iqr = stats::IQR(values)
    )
  }
))

dabg_summaries <- t(apply(
  dabg_matrix,
  2,
  function(values) {
    c(
      dabg_detected_p_lt_0_01 = sum(values < 0.01),
      dabg_detected_p_lt_0_05 = sum(values < 0.05),
      dabg_detected_p_lt_0_10 = sum(values < 0.10),
      dabg_detected_fraction_p_lt_0_01 = mean(values < 0.01),
      dabg_detected_fraction_p_lt_0_05 = mean(values < 0.05),
      dabg_detected_fraction_p_lt_0_10 = mean(values < 0.10),
      dabg_median_p = median(values),
      dabg_mean_p = mean(values)
    )
  }
))

sample_centroid <- colMeans(expression_matrix)
distance_to_centroid <- sqrt(colMeans(
  sweep(expression_matrix, 2, sample_centroid, FUN = "-")^2
))

sample_qc <- data.frame(
  file_name = sample_metadata$file_name,
  animal_id = sample_metadata$animal_id,
  treatment = sample_metadata$treatment,
  age = sample_metadata$age,
  sex = sample_metadata$sex,
  array_prefix = sample_metadata$array_prefix,
  array_position = sample_metadata$array_position,
  PC1 = pca_scores$PC1,
  PC2 = pca_scores$PC2,
  distance_to_centroid = distance_to_centroid,
  expression_summaries,
  dabg_summaries,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_qc$flag_pc1_gt_3sd <- flag_three_sd(abs(sample_qc$PC1), high = TRUE)
sample_qc$flag_pc2_gt_3sd <- flag_three_sd(abs(sample_qc$PC2), high = TRUE)
sample_qc$flag_distance_to_centroid_gt_3sd <- flag_three_sd(
  sample_qc$distance_to_centroid,
  high = TRUE
)
sample_qc$flag_expr_median_gt_3sd <- flag_three_sd(
  sample_qc$expr_median,
  high = TRUE,
  low = TRUE
)
sample_qc$flag_expr_iqr_gt_3sd <- flag_three_sd(
  sample_qc$expr_iqr,
  high = TRUE,
  low = TRUE
)
sample_qc$flag_dabg_fraction_p_lt_0_05_low_3sd <- flag_three_sd(
  sample_qc$dabg_detected_fraction_p_lt_0_05,
  high = FALSE,
  low = TRUE
)

flag_columns <- grep("^flag_", names(sample_qc), value = TRUE)
sample_qc$possible_outlier_flag_count <- rowSums(sample_qc[, flag_columns])
sample_qc$possible_outlier <- sample_qc$possible_outlier_flag_count > 0

write.csv(sample_qc, sample_qc_file, row.names = FALSE, quote = TRUE)

possible_outliers <- sample_qc[sample_qc$possible_outlier, ]
if (nrow(possible_outliers) == 0) {
  possible_outliers <- sample_qc[FALSE, ]
}
write.csv(possible_outliers, outlier_file, row.names = FALSE, quote = TRUE)

summary_lines <- c(
  "miRNA 150001 vendor RMA sample-level QC summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("Expression matrix:", expression_file),
  paste("DABG p-value matrix:", dabg_file),
  paste("Sample sheet:", sample_sheet_file),
  "",
  paste("Samples evaluated:", ncol(expression_matrix)),
  paste("Probesets evaluated:", nrow(expression_matrix)),
  paste("PC1 variance percent:", round(100 * pca_var[1], 2)),
  paste("PC2 variance percent:", round(100 * pca_var[2], 2)),
  paste("PC3 variance percent:", round(100 * pca_var[3], 2)),
  "",
  "Treatment counts:",
  capture.output(print(table(sample_metadata$treatment, useNA = "ifany"))),
  "",
  "Age counts:",
  capture.output(print(table(sample_metadata$age, useNA = "ifany"))),
  "",
  "Sex counts:",
  capture.output(print(table(sample_metadata$sex, useNA = "ifany"))),
  "",
  "Array prefix counts:",
  capture.output(print(table(sample_metadata$array_prefix, useNA = "ifany"))),
  "",
  "Array position counts:",
  capture.output(print(table(sample_metadata$array_position, useNA = "ifany"))),
  "",
  paste("Possible outlier samples flagged:", sum(sample_qc$possible_outlier)),
  "Possible outlier sample IDs:",
  if (any(sample_qc$possible_outlier)) {
    paste(
      sample_qc$animal_id[sample_qc$possible_outlier],
      sample_qc$treatment[sample_qc$possible_outlier],
      paste0("flags=", sample_qc$possible_outlier_flag_count[sample_qc$possible_outlier]),
      sep = " | "
    )
  } else {
    "None"
  },
  "",
  "Outlier flags are review annotations only. No samples were excluded.",
  "Differential expression has not been run."
)

writeLines(summary_lines, summary_file)

message("Wrote annotated PCA scores: ", pca_scores_file)
message("Wrote sample QC metrics: ", sample_qc_file)
message("Wrote possible outlier table: ", outlier_file)
message("Wrote sample QC summary: ", summary_file)
message("Wrote annotated PCA and clustering plots to: ", qc_dir)
