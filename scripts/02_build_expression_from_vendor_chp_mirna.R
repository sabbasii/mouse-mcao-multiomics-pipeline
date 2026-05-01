#!/usr/bin/env Rscript

# Build miRNA 150001 expression matrices from vendor RMA-DABG CHP files.
# This script performs vendor-processed expression extraction and QC only. It
# does not run differential expression.

sample_sheet_file <- file.path(
  "results",
  "sample_sheet",
  "sample_sheet_mirna_150001.csv"
)
chp_dir <- file.path("data", "processed", "mirna_150001", "cc-chp")
qc_root_dir <- file.path("results", "qc")
qc_dir <- file.path(qc_root_dir, "vendor_chp_extraction")
normalized_dir <- file.path("results", "normalized")

expression_file <- file.path(
  normalized_dir,
  "vendor_rma_expression_mirna_150001.csv"
)
dabg_file <- file.path(
  normalized_dir,
  "vendor_dabg_pvalues_mirna_150001.csv"
)
annotation_file <- file.path(
  normalized_dir,
  "vendor_probe_annotation_mirna_150001.csv"
)
matching_summary_file <- file.path(
  qc_dir,
  "vendor_chp_matching_summary_mirna_150001.txt"
)
expression_summary_file <- file.path(
  qc_dir,
  "vendor_chp_expression_summary_mirna_150001.txt"
)

required_sample_columns <- c(
  "platform",
  "file_name",
  "animal_id",
  "sex",
  "age",
  "treatment",
  "is_control",
  "match_status"
)
required_chp_columns <- c(
  "ProbeSetName",
  "QuantificationValue",
  "PValue",
  "ID"
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

make_sample_stem <- function(file_name) {
  sub("[.]CEL$", "", file_name, ignore.case = TRUE)
}

make_chp_stem <- function(file_name) {
  sub("[.]rma-dabg[.]chp$", "", file_name, ignore.case = TRUE)
}

read_vendor_chp <- function(chp_file) {
  chp <- affxparser::readChp(chp_file)
  entries <- chp$QuantificationEntries

  if (is.null(entries)) {
    stop(
      "No QuantificationEntries found in CHP file: ",
      chp_file,
      call. = FALSE
    )
  }

  entries <- as.data.frame(entries, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_chp_columns, names(entries))
  if (length(missing_columns) > 0) {
    stop(
      "CHP file is missing required QuantificationEntries columns: ",
      paste(missing_columns, collapse = ", "),
      "\nFile: ",
      chp_file,
      call. = FALSE
    )
  }

  entries <- entries[, required_chp_columns]
  entries$QuantificationValue <- as.numeric(entries$QuantificationValue)
  entries$PValue <- as.numeric(entries$PValue)
  entries$ID <- as.integer(entries$ID)
  entries
}

summarize_numeric_matrix <- function(matrix_object, label) {
  values <- as.vector(matrix_object)
  finite_values <- values[is.finite(values)]

  lines <- c(
    paste(label, "dimensions:", paste(dim(matrix_object), collapse = " x ")),
    paste(label, "finite values:", sum(is.finite(values))),
    paste(label, "NA values:", sum(is.na(values))),
    paste(label, "NaN values:", sum(is.nan(values))),
    paste(label, "Inf values:", sum(is.infinite(values)))
  )

  if (length(finite_values) > 0) {
    lines <- c(
      lines,
      paste(
        label,
        "finite range:",
        paste(range(finite_values), collapse = " to ")
      ),
      paste(label, "finite summary:"),
      capture.output(print(summary(finite_values)))
    )
  } else {
    lines <- c(lines, paste(label, "finite range: no finite values"))
  }

  lines
}

plot_density_matrix <- function(matrix_object, output_file, main, xlab) {
  finite_columns <- vapply(
    seq_len(ncol(matrix_object)),
    function(index) any(is.finite(matrix_object[, index])),
    logical(1)
  )

  if (!any(finite_columns)) {
    stop("No finite values available for density plot: ", output_file)
  }

  png(output_file, width = 9, height = 7, units = "in", res = 300)
  first_column <- which(finite_columns)[1]
  plot(
    density(matrix_object[, first_column], na.rm = TRUE),
    main = main,
    xlab = xlab,
    col = "gray30",
    lwd = 1
  )
  for (index in which(finite_columns)[-1]) {
    lines(
      density(matrix_object[, index], na.rm = TRUE),
      col = "gray30",
      lwd = 1
    )
  }
  dev.off()
}

require_package("affxparser")

stop_if_missing(sample_sheet_file, "Sample sheet")
stop_if_missing(chp_dir, "CHP directory")

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

chp_files <- list.files(
  chp_dir,
  pattern = "[.]rma-dabg[.]chp$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(chp_files) == 0) {
  stop("No .rma-dabg.chp files found in: ", chp_dir, call. = FALSE)
}

chp_table <- data.frame(
  chp_file = chp_files,
  chp_basename = basename(chp_files),
  chp_stem = make_chp_stem(basename(chp_files)),
  stringsAsFactors = FALSE
)

sample_sheet$sample_stem <- make_sample_stem(sample_sheet$file_name)

if (anyDuplicated(sample_sheet$sample_stem) > 0) {
  duplicated_stems <- unique(
    sample_sheet$sample_stem[duplicated(sample_sheet$sample_stem)]
  )
  stop(
    "Duplicate sample stems in sample sheet: ",
    paste(duplicated_stems, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(chp_table$chp_stem) > 0) {
  duplicated_stems <- unique(chp_table$chp_stem[duplicated(chp_table$chp_stem)])
  stop(
    "Duplicate CHP stems: ",
    paste(duplicated_stems, collapse = ", "),
    call. = FALSE
  )
}

matched_samples <- merge(
  sample_sheet,
  chp_table,
  by.x = "sample_stem",
  by.y = "chp_stem",
  all.x = TRUE,
  sort = FALSE
)

matched_chp <- merge(
  chp_table,
  sample_sheet[, c("sample_stem", "file_name")],
  by.x = "chp_stem",
  by.y = "sample_stem",
  all.x = TRUE,
  sort = FALSE
)

analysis_samples <- matched_samples[
  matched_samples$is_control == FALSE &
    matched_samples$match_status == "matched",
]

if (nrow(analysis_samples) == 0) {
  stop("No matched non-control miRNA samples found.", call. = FALSE)
}

missing_analysis_chp <- analysis_samples$file_name[is.na(analysis_samples$chp_file)]
if (length(missing_analysis_chp) > 0) {
  stop(
    "The following non-control sample-sheet rows have no matching CHP file:\n",
    paste(missing_analysis_chp, collapse = "\n"),
    call. = FALSE
  )
}

analysis_samples <- analysis_samples[order(analysis_samples$file_name), ]
rownames(analysis_samples) <- analysis_samples$file_name

dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(normalized_dir, recursive = TRUE, showWarnings = FALSE)

chp_entries <- lapply(analysis_samples$chp_file, read_vendor_chp)
names(chp_entries) <- analysis_samples$file_name

reference_probe_names <- chp_entries[[1]]$ProbeSetName
reference_probe_ids <- chp_entries[[1]]$ID

for (sample_name in names(chp_entries)) {
  entries <- chp_entries[[sample_name]]

  if (!identical(entries$ProbeSetName, reference_probe_names)) {
    stop(
      "ProbeSetName order/content differs for sample: ",
      sample_name,
      call. = FALSE
    )
  }

  if (!identical(entries$ID, reference_probe_ids)) {
    stop("Probe ID order/content differs for sample: ", sample_name, call. = FALSE)
  }
}

expression_matrix <- do.call(
  cbind,
  lapply(chp_entries, function(entries) entries$QuantificationValue)
)
dabg_matrix <- do.call(
  cbind,
  lapply(chp_entries, function(entries) entries$PValue)
)

rownames(expression_matrix) <- reference_probe_names
rownames(dabg_matrix) <- reference_probe_names
colnames(expression_matrix) <- names(chp_entries)
colnames(dabg_matrix) <- names(chp_entries)

probe_annotation <- data.frame(
  ProbeSetName = reference_probe_names,
  ID = reference_probe_ids,
  stringsAsFactors = FALSE
)

write.csv(
  data.frame(
    ProbeSetName = rownames(expression_matrix),
    expression_matrix,
    check.names = FALSE
  ),
  expression_file,
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  data.frame(
    ProbeSetName = rownames(dabg_matrix),
    dabg_matrix,
    check.names = FALSE
  ),
  dabg_file,
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  probe_annotation,
  annotation_file,
  row.names = FALSE,
  quote = TRUE
)

png(
  file.path(qc_dir, "vendor_rma_boxplot_mirna_150001.png"),
  width = 11,
  height = 7,
  units = "in",
  res = 300
)
boxplot(
  expression_matrix,
  las = 2,
  main = "miRNA 150001 vendor RMA-DABG expression",
  ylab = "Vendor RMA log2 expression",
  cex.axis = 0.45
)
dev.off()

plot_density_matrix(
  expression_matrix,
  file.path(qc_dir, "vendor_rma_density_mirna_150001.png"),
  "miRNA 150001 vendor RMA-DABG density",
  "Vendor RMA log2 expression"
)

complete_expression_rows <- stats::complete.cases(expression_matrix) &
  apply(expression_matrix, 1, function(values) all(is.finite(values)))

if (sum(complete_expression_rows) < 2) {
  stop("Not enough complete finite probesets for PCA/clustering.", call. = FALSE)
}

pca <- prcomp(t(expression_matrix[complete_expression_rows, ]), center = TRUE)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)

png(
  file.path(qc_dir, "vendor_rma_pca_mirna_150001.png"),
  width = 8,
  height = 7,
  units = "in",
  res = 300
)
plot(
  pca$x[, 1],
  pca$x[, 2],
  pch = 19,
  col = as.integer(factor(analysis_samples$treatment)),
  xlab = paste0("PC1 (", round(100 * pca_var[1], 1), "%)"),
  ylab = paste0("PC2 (", round(100 * pca_var[2], 1), "%)"),
  main = "miRNA 150001 PCA from vendor RMA-DABG"
)
legend(
  "topright",
  legend = levels(factor(analysis_samples$treatment)),
  col = seq_along(levels(factor(analysis_samples$treatment))),
  pch = 19,
  cex = 0.8
)
dev.off()

sample_dist <- dist(t(expression_matrix[complete_expression_rows, ]))
sample_hclust <- hclust(sample_dist, method = "complete")

png(
  file.path(qc_dir, "vendor_rma_sample_clustering_mirna_150001.png"),
  width = 11,
  height = 7,
  units = "in",
  res = 300
)
plot(
  sample_hclust,
  labels = analysis_samples$animal_id,
  main = "miRNA 150001 sample clustering from vendor RMA-DABG",
  xlab = "",
  sub = ""
)
dev.off()

matching_summary <- c(
  "miRNA 150001 vendor CHP sample matching summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("Sample sheet:", sample_sheet_file),
  paste("Input CHP directory:", chp_dir),
  paste("RMA-DABG CHP files found:", nrow(chp_table)),
  paste("Sample sheet rows:", nrow(sample_sheet)),
  paste(
    "Sample sheet rows with matching CHP:",
    sum(!is.na(matched_samples$chp_file))
  ),
  paste(
    "Sample sheet rows without matching CHP:",
    sum(is.na(matched_samples$chp_file))
  ),
  paste(
    "CHP files without sample-sheet row:",
    sum(is.na(matched_chp$file_name))
  ),
  paste("Control rows excluded from main matrices:", sum(sample_sheet$is_control)),
  paste("Matched non-control samples used:", nrow(analysis_samples)),
  "",
  "Unmatched sample-sheet rows:",
  if (any(is.na(matched_samples$chp_file))) {
    matched_samples$file_name[is.na(matched_samples$chp_file)]
  } else {
    "None"
  },
  "",
  "Unmatched CHP files:",
  if (any(is.na(matched_chp$file_name))) {
    matched_chp$chp_basename[is.na(matched_chp$file_name)]
  } else {
    "None"
  }
)

expression_summary <- c(
  "miRNA 150001 vendor CHP expression summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("Expression matrix output:", expression_file),
  paste("DABG p-value matrix output:", dabg_file),
  paste("Probe annotation output:", annotation_file),
  "",
  summarize_numeric_matrix(expression_matrix, "Expression matrix"),
  "",
  summarize_numeric_matrix(dabg_matrix, "DABG p-value matrix"),
  "",
  paste("Probesets:", nrow(expression_matrix)),
  paste("Samples:", ncol(expression_matrix)),
  paste("Complete finite probesets used for PCA/clustering:", sum(complete_expression_rows)),
  "",
  "Treatment counts used:",
  capture.output(print(table(analysis_samples$treatment, useNA = "ifany"))),
  "",
  "Age counts used:",
  capture.output(print(table(analysis_samples$age, useNA = "ifany"))),
  "",
  "Sex counts used:",
  capture.output(print(table(analysis_samples$sex, useNA = "ifany")))
)

writeLines(matching_summary, matching_summary_file)
writeLines(expression_summary, expression_summary_file)

message("Wrote vendor expression matrix: ", expression_file)
message("Wrote vendor DABG p-value matrix: ", dabg_file)
message("Wrote vendor probe annotation: ", annotation_file)
message("Wrote vendor CHP matching summary: ", matching_summary_file)
message("Wrote vendor expression summary: ", expression_summary_file)
message("Wrote vendor QC plots to: ", qc_dir)
