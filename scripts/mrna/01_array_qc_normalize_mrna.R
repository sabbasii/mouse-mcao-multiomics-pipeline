#!/usr/bin/env Rscript

# Read, inspect, and RMA-normalize Clariom S/mRNA 150002 CEL files.
# This script performs array-level preprocessing only. It does not run
# differential expression or biological interpretation.
#
# Prerequisite:
#   Rscript scripts/mrna/00_build_sample_sheet_mrna.R
#
# Run from the repository root:
#   Rscript scripts/mrna/01_array_qc_normalize_mrna.R

sample_sheet_file <- file.path(
  "results",
  "mrna",
  "sample_sheet",
  "sample_sheet_mrna_150002.csv"
)

qc_dir <- file.path("results", "mrna", "qc", "cel_rma")
normalized_dir <- file.path("results", "mrna", "normalized")

raw_boxplot_file <- file.path(
  qc_dir,
  "raw_cel_intensity_boxplot_mrna_150002.png"
)
raw_density_file <- file.path(
  qc_dir,
  "raw_cel_intensity_density_mrna_150002.png"
)
normalized_boxplot_file <- file.path(
  qc_dir,
  "rma_expression_boxplot_mrna_150002.png"
)
normalized_density_file <- file.path(
  qc_dir,
  "rma_expression_density_mrna_150002.png"
)
qc_metrics_file <- file.path(
  qc_dir,
  "cel_rma_sample_metrics_mrna_150002.csv"
)
summary_file <- file.path(
  qc_dir,
  "cel_rma_summary_mrna_150002.txt"
)

rma_matrix_file <- file.path(
  normalized_dir,
  "rma_expression_mrna_150002.csv"
)
rma_eset_file <- file.path(
  normalized_dir,
  "rma_expression_set_mrna_150002.rds"
)

required_packages <- c("Biobase", "limma", "oligo")
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them in the project renv environment before continuing.",
    call. = FALSE
  )
}

if (!file.exists(sample_sheet_file)) {
  stop(
    "Sample sheet not found: ",
    sample_sheet_file,
    ". Run scripts/mrna/00_build_sample_sheet_mrna.R first.",
    call. = FALSE
  )
}

sample_sheet <- read.csv(
  sample_sheet_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_sample_columns <- c(
  "file_name",
  "cel_path",
  "animal_id",
  "is_control",
  "match_status"
)
missing_sample_columns <- setdiff(
  required_sample_columns,
  names(sample_sheet)
)

if (length(missing_sample_columns) > 0L) {
  stop(
    "Sample sheet is missing required columns: ",
    paste(missing_sample_columns, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(sample_sheet$file_name)) {
  stop("Sample sheet contains duplicate CEL filenames.", call. = FALSE)
}

allowed_statuses <- c("matched", "control_without_animal_metadata")
unexpected_statuses <- setdiff(
  unique(sample_sheet$match_status),
  allowed_statuses
)

if (length(unexpected_statuses) > 0L) {
  stop(
    "Sample sheet contains unresolved match status(es): ",
    paste(unexpected_statuses, collapse = ", "),
    call. = FALSE
  )
}

missing_cel_files <- sample_sheet$cel_path[
  !file.exists(sample_sheet$cel_path)
]
if (length(missing_cel_files) > 0L) {
  stop(
    "CEL file(s) listed in the sample sheet were not found: ",
    paste(missing_cel_files, collapse = ", "),
    call. = FALSE
  )
}

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(normalized_dir, recursive = TRUE, showWarnings = FALSE)

message("Reading ", nrow(sample_sheet), " Clariom S CEL files...")
raw_cel <- oligo::read.celfiles(
  sample_sheet$cel_path,
  verbose = TRUE
)

raw_intensities <- Biobase::exprs(raw_cel)

if (!all(is.finite(raw_intensities))) {
  stop(
    "Raw CEL intensities contain NA, NaN, or infinite values.",
    call. = FALSE
  )
}

cel_object_names <- basename(Biobase::sampleNames(raw_cel))
if (!identical(cel_object_names, sample_sheet$file_name)) {
  stop(
    "CEL order in the loaded object does not match the sample sheet.",
    call. = FALSE
  )
}

platform_package <- Biobase::annotation(raw_cel)
message("Detected platform-design package: ", platform_package)

png(
  raw_boxplot_file,
  width = 2400,
  height = 1400,
  res = 180
)
boxplot(
  log2(raw_intensities),
  outline = FALSE,
  las = 2,
  names = sub("\\.[Cc][Ee][Ll]$", "", sample_sheet$file_name),
  main = "Raw Clariom S CEL intensity distributions",
  ylab = "log2 raw intensity"
)
dev.off()

png(
  raw_density_file,
  width = 2200,
  height = 1400,
  res = 180
)
limma::plotDensities(
  log2(raw_intensities),
  legend = FALSE,
  main = "Raw Clariom S CEL intensity densities",
  xlab = "log2 raw intensity"
)
dev.off()

raw_metrics <- data.frame(
  file_name = sample_sheet$file_name,
  raw_min = apply(raw_intensities, 2L, min),
  raw_q1 = apply(
    raw_intensities,
    2L,
    stats::quantile,
    probs = 0.25,
    names = FALSE
  ),
  raw_median = apply(raw_intensities, 2L, stats::median),
  raw_mean = colMeans(raw_intensities),
  raw_q3 = apply(
    raw_intensities,
    2L,
    stats::quantile,
    probs = 0.75,
    names = FALSE
  ),
  raw_max = apply(raw_intensities, 2L, max),
  stringsAsFactors = FALSE
)

message("Running oligo::rma()...")
rma_eset <- oligo::rma(raw_cel)
rma_matrix <- Biobase::exprs(rma_eset)

rma_na_count <- sum(is.na(rma_matrix))
rma_nan_count <- sum(is.nan(rma_matrix))
rma_inf_count <- sum(is.infinite(rma_matrix))
rma_finite_count <- sum(is.finite(rma_matrix))

if (!all(is.finite(rma_matrix))) {
  stop(
    "RMA produced invalid expression values. ",
    "Finite: ", rma_finite_count,
    "; NA: ", rma_na_count,
    "; NaN: ", rma_nan_count,
    "; Inf: ", rma_inf_count,
    ". No normalized outputs were written.",
    call. = FALSE
  )
}

if (!identical(
  basename(colnames(rma_matrix)),
  sample_sheet$file_name
)) {
  stop(
    "RMA expression columns do not match the sample-sheet order.",
    call. = FALSE
  )
}

png(
  normalized_boxplot_file,
  width = 2400,
  height = 1400,
  res = 180
)
boxplot(
  rma_matrix,
  outline = FALSE,
  las = 2,
  names = sub("\\.[Cc][Ee][Ll]$", "", sample_sheet$file_name),
  main = "RMA-normalized Clariom S expression distributions",
  ylab = "RMA expression (log2)"
)
dev.off()

png(
  normalized_density_file,
  width = 2200,
  height = 1400,
  res = 180
)
limma::plotDensities(
  rma_matrix,
  legend = FALSE,
  main = "RMA-normalized Clariom S expression densities",
  xlab = "RMA expression (log2)"
)
dev.off()

normalized_metrics <- data.frame(
  file_name = sample_sheet$file_name,
  rma_min = apply(rma_matrix, 2L, min),
  rma_q1 = apply(
    rma_matrix,
    2L,
    stats::quantile,
    probs = 0.25,
    names = FALSE
  ),
  rma_median = apply(rma_matrix, 2L, stats::median),
  rma_mean = colMeans(rma_matrix),
  rma_q3 = apply(
    rma_matrix,
    2L,
    stats::quantile,
    probs = 0.75,
    names = FALSE
  ),
  rma_max = apply(rma_matrix, 2L, max),
  stringsAsFactors = FALSE
)

qc_metrics <- merge(
  sample_sheet,
  raw_metrics,
  by = "file_name",
  all.x = TRUE,
  sort = FALSE
)
qc_metrics <- merge(
  qc_metrics,
  normalized_metrics,
  by = "file_name",
  all.x = TRUE,
  sort = FALSE
)
qc_metrics <- qc_metrics[
  match(sample_sheet$file_name, qc_metrics$file_name),
  ,
  drop = FALSE
]

rma_output <- data.frame(
  transcript_cluster_id = rownames(rma_matrix),
  rma_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

write.csv(
  qc_metrics,
  qc_metrics_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  rma_output,
  rma_matrix_file,
  row.names = FALSE,
  na = ""
)
saveRDS(rma_eset, rma_eset_file)

summary_lines <- c(
  "Clariom S/mRNA 150002 CEL and RMA summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("Sample sheet:", sample_sheet_file),
  paste("Arrays read:", ncol(raw_intensities)),
  paste("Raw probe features per array:", nrow(raw_intensities)),
  paste("Raw finite values:", sum(is.finite(raw_intensities))),
  paste("Raw NA values:", sum(is.na(raw_intensities))),
  paste("Raw NaN values:", sum(is.nan(raw_intensities))),
  paste("Raw infinite values:", sum(is.infinite(raw_intensities))),
  paste("Detected platform-design package:", platform_package),
  "",
  paste("RMA transcript clusters:", nrow(rma_matrix)),
  paste("RMA arrays:", ncol(rma_matrix)),
  paste("RMA finite values:", rma_finite_count),
  paste("RMA NA values:", rma_na_count),
  paste("RMA NaN values:", rma_nan_count),
  paste("RMA infinite values:", rma_inf_count),
  "",
  paste("RMA expression matrix:", rma_matrix_file),
  paste("RMA ExpressionSet:", rma_eset_file),
  paste("Sample QC metrics:", qc_metrics_file),
  paste("Raw boxplot:", raw_boxplot_file),
  paste("Raw density plot:", raw_density_file),
  paste("Normalized boxplot:", normalized_boxplot_file),
  paste("Normalized density plot:", normalized_density_file)
)

writeLines(summary_lines, summary_file)

message("Wrote RMA expression matrix: ", rma_matrix_file)
message("Wrote RMA ExpressionSet: ", rma_eset_file)
message("Wrote CEL/RMA QC outputs under: ", qc_dir)
