#!/usr/bin/env Rscript

# Array-level QC and normalization for miRNA 150001.
# This script does not run differential expression.

sample_sheet_file <- file.path(
  "results",
  "sample_sheet",
  "sample_sheet_mirna_150001.csv"
)
qc_dir <- file.path("results", "qc")
normalized_dir <- file.path("results", "normalized")

required_columns <- c(
  "platform",
  "file_name",
  "cel_path",
  "animal_id",
  "sex",
  "age",
  "treatment",
  "is_control",
  "match_status"
)

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Required R package is not installed: ",
      package,
      "\nInstall it before running this script.",
      call. = FALSE
    )
  }
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

require_package("oligo")
require_package("Biobase")

stop_if_missing(sample_sheet_file, "Sample sheet")

sample_sheet <- read.csv(
  sample_sheet_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_columns <- setdiff(required_columns, names(sample_sheet))
if (length(missing_columns) > 0) {
  stop(
    "Sample sheet is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

analysis_samples <- sample_sheet[
  sample_sheet$is_control == FALSE &
    sample_sheet$match_status == "matched",
]

if (nrow(analysis_samples) == 0) {
  stop("No matched non-control miRNA samples found.", call. = FALSE)
}

missing_cel <- analysis_samples$cel_path[!file.exists(analysis_samples$cel_path)]
if (length(missing_cel) > 0) {
  stop(
    "The following CEL files are listed in the sample sheet but missing:\n",
    paste(missing_cel, collapse = "\n"),
    call. = FALSE
  )
}

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(normalized_dir, recursive = TRUE, showWarnings = FALSE)

cel_paths <- analysis_samples$cel_path
names(cel_paths) <- analysis_samples$file_name

raw_data <- oligo::read.celfiles(cel_paths)
sample_names <- analysis_samples$file_name
rownames(analysis_samples) <- sample_names
Biobase::sampleNames(raw_data) <- sample_names
Biobase::pData(raw_data) <- analysis_samples
rownames(Biobase::pData(Biobase::protocolData(raw_data))) <- sample_names

platform_info <- list(
  annotation = Biobase::annotation(raw_data),
  manufacturer = tryCatch(
    unique(Biobase::protocolData(raw_data)$manufacturer),
    error = function(error) character()
  ),
  assay_type = tryCatch(
    unique(Biobase::protocolData(raw_data)$assayType),
    error = function(error) character()
  ),
  array_type = tryCatch(
    unique(Biobase::protocolData(raw_data)$arrayType),
    error = function(error) character()
  )
)

platform_lines <- c(
  "miRNA 150001 platform check",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("CEL files inspected:", length(cel_paths)),
  paste("annotation:", paste(platform_info$annotation, collapse = "; ")),
  paste("manufacturer:", paste(platform_info$manufacturer, collapse = "; ")),
  paste("assay_type:", paste(platform_info$assay_type, collapse = "; ")),
  paste("array_type:", paste(platform_info$array_type, collapse = "; "))
)

writeLines(platform_lines, file.path(qc_dir, "mirna_150001_platform_check.txt"))

identified_platform <- any(
  nzchar(unlist(platform_info, use.names = FALSE))
)

if (!identified_platform) {
  stop(
    "Unable to identify CEL platform metadata. ",
    "Wrote platform check to ",
    file.path(qc_dir, "mirna_150001_platform_check.txt"),
    ". Normalization was not run.",
    call. = FALSE
  )
}

normalized_data <- oligo::rma(raw_data)
expr <- Biobase::exprs(normalized_data)

Biobase::sampleNames(normalized_data) <- sample_names
colnames(expr) <- sample_names
Biobase::pData(normalized_data) <- analysis_samples
rownames(Biobase::pData(Biobase::protocolData(normalized_data))) <- sample_names

saveRDS(
  normalized_data,
  file.path(normalized_dir, "normalized_expression_mirna_150001_eset.rds")
)

write.csv(
  expr,
  file.path(normalized_dir, "normalized_expression_mirna_150001.csv"),
  quote = TRUE
)

write.csv(
  analysis_samples,
  file.path(normalized_dir, "analysis_samples_mirna_150001.csv"),
  row.names = FALSE,
  quote = TRUE
)

pdf(file.path(qc_dir, "mirna_150001_raw_boxplot.pdf"), width = 11, height = 7)
boxplot(
  raw_data,
  target = "core",
  las = 2,
  main = "miRNA 150001 raw CEL intensities",
  cex.axis = 0.45
)
dev.off()

pdf(file.path(qc_dir, "mirna_150001_normalized_boxplot.pdf"), width = 11, height = 7)
boxplot(
  expr,
  las = 2,
  main = "miRNA 150001 RMA-normalized expression",
  ylab = "log2 expression",
  cex.axis = 0.45
)
dev.off()

pdf(file.path(qc_dir, "mirna_150001_normalized_density.pdf"), width = 9, height = 7)
plot(
  density(expr[, 1]),
  main = "miRNA 150001 RMA-normalized density",
  xlab = "log2 expression",
  col = "gray30",
  lwd = 1
)
if (ncol(expr) > 1) {
  for (i in seq_len(ncol(expr))[-1]) {
    lines(density(expr[, i]), col = "gray30", lwd = 1)
  }
}
dev.off()

pca <- prcomp(t(expr), center = TRUE, scale. = FALSE)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
pca_plot <- data.frame(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  sex = analysis_samples$sex,
  age = analysis_samples$age,
  treatment = analysis_samples$treatment,
  stringsAsFactors = FALSE
)

write.csv(
  pca_plot,
  file.path(qc_dir, "mirna_150001_normalized_pca_scores.csv"),
  row.names = FALSE,
  quote = TRUE
)

pdf(file.path(qc_dir, "mirna_150001_normalized_pca.pdf"), width = 8, height = 7)
plot(
  pca_plot$PC1,
  pca_plot$PC2,
  pch = 19,
  col = as.integer(factor(pca_plot$treatment)),
  xlab = paste0("PC1 (", round(100 * pca_var[1], 1), "%)"),
  ylab = paste0("PC2 (", round(100 * pca_var[2], 1), "%)"),
  main = "miRNA 150001 PCA after RMA normalization"
)
legend(
  "topright",
  legend = levels(factor(pca_plot$treatment)),
  col = seq_along(levels(factor(pca_plot$treatment))),
  pch = 19,
  cex = 0.8
)
dev.off()

sample_dist <- dist(t(expr))
sample_hclust <- hclust(sample_dist, method = "complete")

pdf(file.path(qc_dir, "mirna_150001_normalized_clustering.pdf"), width = 11, height = 7)
plot(
  sample_hclust,
  labels = analysis_samples$animal_id,
  main = "miRNA 150001 sample clustering after RMA normalization",
  xlab = "",
  sub = ""
)
dev.off()

qc_summary <- c(
  "miRNA 150001 QC and normalization summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("Sample sheet:", sample_sheet_file),
  paste("Input rows:", nrow(sample_sheet)),
  paste("Matched non-control samples used:", nrow(analysis_samples)),
  paste("Controls excluded:", sum(sample_sheet$is_control == TRUE, na.rm = TRUE)),
  paste(
    "Non-matched rows excluded:",
    sum(sample_sheet$match_status != "matched", na.rm = TRUE)
  ),
  "",
  "Treatment counts used:",
  capture.output(print(table(analysis_samples$treatment, useNA = "ifany"))),
  "",
  "Age counts used:",
  capture.output(print(table(analysis_samples$age, useNA = "ifany"))),
  "",
  "Sex counts used:",
  capture.output(print(table(analysis_samples$sex, useNA = "ifany"))),
  "",
  paste(
    "Normalized expression CSV:",
    file.path(normalized_dir, "normalized_expression_mirna_150001.csv")
  ),
  paste(
    "Normalized ExpressionSet RDS:",
    file.path(normalized_dir, "normalized_expression_mirna_150001_eset.rds")
  )
)

writeLines(qc_summary, file.path(qc_dir, "mirna_150001_qc_normalization_summary.txt"))

message("Wrote QC outputs to: ", qc_dir)
message("Wrote normalized outputs to: ", normalized_dir)
