#!/usr/bin/env Rscript

# Sample-level QC for one prepared miRNA expression dataset.
# Usage: Rscript scripts/mirna/03_sample_qc_mirna.R <dataset>

args <- commandArgs(trailingOnly = TRUE)
valid_datasets <- c(
  "vendor_chp", "complete_rma_excluding_two", "rma_normalized_mirna"
)
if (length(args) != 1L || !args[[1]] %in% valid_datasets) {
  stop(
    "Supply exactly one dataset: ",
    paste(valid_datasets, collapse = ", "),
    call. = FALSE
  )
}
dataset <- args[[1]]

expression_file <- file.path(
  "results", "mirna", "expression", dataset,
  "expression_matrix_mirna.csv"
)
sample_sheet_file <- file.path(
  "results", "mirna", "sample_sheet", "sample_sheet_mirna.csv"
)
dabg_file <- file.path(
  "results", "mirna", "expression", "vendor_chp",
  "dabg_pvalues_mirna.csv"
)
output_dir <- file.path(
  "results", "mirna", "qc", dataset, "sample_qc"
)

read_probe_matrix <- function(path) {
  table <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"ProbeSetName" %in% names(table)) {
    stop("Missing ProbeSetName column: ", path, call. = FALSE)
  }
  sample_columns <- setdiff(
    names(table), c("ProbeSetName", "feature_name", "ID")
  )
  matrix_object <- as.matrix(table[, sample_columns, drop = FALSE])
  storage.mode(matrix_object) <- "numeric"
  rownames(matrix_object) <- as.character(table$ProbeSetName)
  colnames(matrix_object) <- sub(
    "[.]CEL$", "", colnames(matrix_object), ignore.case = TRUE
  )
  matrix_object
}

flag_three_sd <- function(values, direction = c("high", "low")) {
  direction <- match.arg(direction)
  center <- mean(values, na.rm = TRUE)
  spread <- sd(values, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) {
    return(rep(FALSE, length(values)))
  }
  if (direction == "high") {
    values > center + 3 * spread
  } else {
    values < center - 3 * spread
  }
}

plot_pca <- function(scores, variable, path) {
  group <- factor(scores[[variable]])
  colors <- grDevices::hcl.colors(max(3L, nlevels(group)), "Dark 3")[
    as.integer(group)
  ]
  png(path, width = 9, height = 7, units = "in", res = 200)
  plot(
    scores$PC1, scores$PC2,
    col = colors, pch = 19,
    xlab = paste0("PC1 (", scores$PC1_percent[[1]], "%)"),
    ylab = paste0("PC2 (", scores$PC2_percent[[1]], "%)"),
    main = paste(dataset, "PCA by", variable)
  )
  text(scores$PC1, scores$PC2, labels = scores$animal_id, pos = 3, cex = 0.55)
  legend(
    "topright", legend = levels(group),
    col = grDevices::hcl.colors(max(3L, nlevels(group)), "Dark 3")[
      seq_len(nlevels(group))
    ],
    pch = 19, cex = 0.8
  )
  dev.off()
}

if (!file.exists(expression_file) || !file.exists(sample_sheet_file)) {
  stop("Required expression matrix or sample sheet is missing.", call. = FALSE)
}

expression <- read_probe_matrix(expression_file)
input_dimensions <- dim(expression)
available_samples <- colnames(expression)
missing_count <- sum(is.na(expression))
nan_count <- sum(is.nan(expression))
infinite_count <- sum(is.infinite(expression))
if (any(!is.finite(expression))) {
  stop("Expression matrix contains non-finite values.", call. = FALSE)
}
sample_sheet <- read.csv(
  sample_sheet_file, stringsAsFactors = FALSE, check.names = FALSE
)
sample_sheet$sample_name <- sub(
  "[.]CEL$", "", sample_sheet$file_name, ignore.case = TRUE
)
eligible <- sample_sheet$match_status == "matched" & !sample_sheet$is_control
metadata <- sample_sheet[eligible, , drop = FALSE]
metadata <- metadata[metadata$sample_name %in% colnames(expression), , drop = FALSE]
if (nrow(metadata) < 3L) {
  stop("Fewer than three biological samples overlap expression and metadata.")
}
expression <- expression[, metadata$sample_name, drop = FALSE]

pca <- prcomp(t(expression), center = TRUE, scale. = FALSE)
variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
scores <- cbind(
  metadata,
  data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC1_percent = variance[[1]],
    PC2_percent = variance[[2]],
    row.names = NULL
  )
)
scores$array_prefix <- sub("_.*$", "", scores$sample_name)
scores$pca_distance <- sqrt(scores$PC1^2 + scores$PC2^2)

sample_metrics <- data.frame(
  sample_name = colnames(expression),
  expression_mean = colMeans(expression),
  expression_median = apply(expression, 2, median),
  expression_sd = apply(expression, 2, sd),
  distance_from_centroid = sqrt(
    colSums((expression - rowMeans(expression))^2)
  ),
  stringsAsFactors = FALSE
)
sample_metrics <- merge(
  metadata, sample_metrics, by = "sample_name", sort = FALSE
)
sample_metrics$pca_distance <- scores$pca_distance[
  match(sample_metrics$sample_name, scores$sample_name)
]
sample_metrics$flag_expression_distance <- flag_three_sd(
  sample_metrics$distance_from_centroid, "high"
)
sample_metrics$flag_pca_distance <- flag_three_sd(
  sample_metrics$pca_distance, "high"
)

if (dataset == "vendor_chp" && file.exists(dabg_file)) {
  dabg <- read_probe_matrix(dabg_file)
  dabg <- dabg[, sample_metrics$sample_name, drop = FALSE]
  sample_metrics$dabg_fraction_p_lt_0_05 <- colMeans(dabg < 0.05)
  sample_metrics$flag_low_dabg_detection <- flag_three_sd(
    sample_metrics$dabg_fraction_p_lt_0_05, "low"
  )
}
flag_columns <- grep("^flag_", names(sample_metrics), value = TRUE)
sample_metrics$possible_outlier_review <- rowSums(
  sample_metrics[, flag_columns, drop = FALSE]
) > 0

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

png(
  file.path(output_dir, "expression_boxplot_mirna.png"),
  width = 14, height = 8, units = "in", res = 200
)
boxplot(
  expression,
  names = metadata$animal_id,
  las = 2,
  outline = FALSE,
  cex.axis = 0.65,
  main = paste(dataset, "expression distributions"),
  ylab = "Expression (log2)"
)
dev.off()

density_colors <- grDevices::adjustcolor(
  grDevices::hcl.colors(ncol(expression), "Dynamic"),
  alpha.f = 0.65
)
png(
  file.path(output_dir, "expression_density_mirna.png"),
  width = 10, height = 7, units = "in", res = 200
)
first_density <- density(expression[, 1])
plot(
  first_density,
  col = density_colors[[1]],
  lwd = 1,
  main = paste(dataset, "expression densities"),
  xlab = "Expression (log2)"
)
if (ncol(expression) > 1L) {
  for (index in 2:ncol(expression)) {
    lines(
      density(expression[, index]),
      col = density_colors[[index]],
      lwd = 1
    )
  }
}
dev.off()

write.csv(
  scores,
  file.path(output_dir, "pca_scores_mirna.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  sample_metrics,
  file.path(output_dir, "sample_qc_metrics_mirna.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  sample_metrics[sample_metrics$possible_outlier_review, , drop = FALSE],
  file.path(output_dir, "possible_outliers_mirna.csv"),
  row.names = FALSE,
  quote = TRUE
)

for (variable in c("treatment", "age", "sex", "array_prefix")) {
  plot_pca(
    scores, variable,
    file.path(output_dir, paste0("pca_by_", variable, "_mirna.png"))
  )
}

png(
  file.path(output_dir, "sample_clustering_mirna.png"),
  width = 11, height = 9, units = "in", res = 200
)
plot(
  hclust(dist(t(expression))),
  labels = metadata$animal_id,
  main = paste(dataset, "sample clustering"),
  xlab = "", sub = ""
)
dev.off()

writeLines(
  c(
    paste("Dataset:", dataset),
    paste(
      "Input expression dimensions:",
      paste(input_dimensions, collapse = " x ")
    ),
    paste(
      "QC expression dimensions:",
      paste(dim(expression), collapse = " x ")
    ),
    paste("Missing values:", missing_count),
    paste("NaN values:", nan_count),
    paste("Infinite values:", infinite_count),
    paste("Probesets used for QC:", nrow(expression), "(all available probesets)"),
    paste("Biological samples used:", ncol(expression)),
    paste("Metadata samples absent from this dataset:", sum(
      eligible & !sample_sheet$sample_name %in% available_samples
    )),
    paste("PC1 variance (%):", round(variance[[1]], 2)),
    paste("PC2 variance (%):", round(variance[[2]], 2)),
    paste(
      "Samples flagged for review:",
      sum(sample_metrics$possible_outlier_review)
    ),
    "No samples are excluded automatically by this QC script."
  ),
  file.path(output_dir, "sample_qc_summary_mirna.txt")
)

message("Completed sample QC for: ", dataset)
message("QC output: ", output_dir)
