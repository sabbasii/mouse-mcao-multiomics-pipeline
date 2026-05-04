#!/usr/bin/env Rscript

# First-pass miRNA 150001 differential expression with limma.
# This treatment-only analysis is exploratory. It does not include covariates,
# pathway analysis, or sample exclusion based on QC flags.

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
de_root_dir <- file.path("results", "differential_expression")
output_dir <- file.path(de_root_dir, "first_pass_limma")

included_treatments <- c("Sham", "MCAO1hr", "MCAO3hr")
contrast_definitions <- c(
  MCAO1hr_vs_Sham = "MCAO1hr - Sham",
  MCAO3hr_vs_Sham = "MCAO3hr - Sham"
)
dabg_detection_p <- 0.05
dabg_detection_fraction <- 0.20

sample_table_file <- file.path(
  output_dir,
  "first_pass_sample_table_mirna_150001.csv"
)
design_file <- file.path(
  output_dir,
  "first_pass_design_matrix_mirna_150001.csv"
)
contrast_file <- file.path(
  output_dir,
  "first_pass_contrast_matrix_mirna_150001.csv"
)
dabg_summary_file <- file.path(
  output_dir,
  "first_pass_dabg_detection_summary_mirna_150001.csv"
)
summary_file <- file.path(
  output_dir,
  "first_pass_limma_summary_mirna_150001.txt"
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

write_limma_results <- function(fit, coefficient, output_file, detection_summary) {
  result <- limma::topTable(
    fit,
    coef = coefficient,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )

  result$ProbeSetName <- rownames(result)
  result <- merge(
    result,
    detection_summary,
    by = "ProbeSetName",
    all.x = TRUE,
    sort = FALSE
  )

  result <- result[
    ,
    c(
      "ProbeSetName",
      setdiff(names(result), "ProbeSetName")
    )
  ]

  write.csv(result, output_file, row.names = FALSE, quote = TRUE)
  result
}

require_package("limma")

stop_if_missing(expression_file, "Vendor RMA expression matrix")
stop_if_missing(dabg_file, "Vendor DABG p-value matrix")
stop_if_missing(sample_sheet_file, "Sample sheet")

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

analysis_samples <- sample_sheet[
  sample_sheet$is_control == FALSE &
    sample_sheet$match_status == "matched" &
    sample_sheet$treatment %in% included_treatments,
]

if (nrow(analysis_samples) == 0) {
  stop("No samples remain after first-pass inclusion rules.", call. = FALSE)
}

missing_expression_samples <- setdiff(
  analysis_samples$file_name,
  colnames(expression_matrix)
)
if (length(missing_expression_samples) > 0) {
  stop(
    "Included samples are missing from expression matrix:\n",
    paste(missing_expression_samples, collapse = "\n"),
    call. = FALSE
  )
}

expression_matrix <- expression_matrix[, analysis_samples$file_name, drop = FALSE]
dabg_matrix <- dabg_matrix[, analysis_samples$file_name, drop = FALSE]

metadata_match <- match(colnames(expression_matrix), analysis_samples$file_name)
if (any(is.na(metadata_match))) {
  stop(
    "Expression matrix columns have no matching included metadata row:\n",
    paste(colnames(expression_matrix)[is.na(metadata_match)], collapse = "\n"),
    call. = FALSE
  )
}

analysis_samples <- analysis_samples[metadata_match, ]
analysis_samples <- analysis_samples[order(analysis_samples$treatment), ]
rownames(analysis_samples) <- analysis_samples$file_name

expression_matrix <- expression_matrix[, analysis_samples$file_name, drop = FALSE]
dabg_matrix <- dabg_matrix[, analysis_samples$file_name, drop = FALSE]

treatment <- factor(
  analysis_samples$treatment,
  levels = c("Sham", "MCAO1hr", "MCAO3hr")
)
design <- model.matrix(~ 0 + treatment)
colnames(design) <- levels(treatment)
rownames(design) <- analysis_samples$file_name

contrast_matrix <- limma::makeContrasts(
  contrasts = contrast_definitions,
  levels = design
)

dabg_detected_p_lt_0_01 <- rowSums(dabg_matrix < 0.01)
dabg_detected_p_lt_0_05 <- rowSums(dabg_matrix < dabg_detection_p)
dabg_detected_p_lt_0_10 <- rowSums(dabg_matrix < 0.10)
dabg_detection_summary <- data.frame(
  ProbeSetName = rownames(dabg_matrix),
  dabg_detected_p_lt_0_01 = dabg_detected_p_lt_0_01,
  dabg_detected_p_lt_0_05 = dabg_detected_p_lt_0_05,
  dabg_detected_p_lt_0_10 = dabg_detected_p_lt_0_10,
  dabg_detected_fraction_p_lt_0_01 = dabg_detected_p_lt_0_01 / ncol(dabg_matrix),
  dabg_detected_fraction_p_lt_0_05 = dabg_detected_p_lt_0_05 / ncol(dabg_matrix),
  dabg_detected_fraction_p_lt_0_10 = dabg_detected_p_lt_0_10 / ncol(dabg_matrix),
  dabg_light_filter_keep = dabg_detected_p_lt_0_05 >=
    ceiling(dabg_detection_fraction * ncol(dabg_matrix)),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(analysis_samples, sample_table_file, row.names = FALSE, quote = TRUE)
write.csv(design, design_file, quote = TRUE)
write.csv(contrast_matrix, contrast_file, quote = TRUE)
write.csv(dabg_detection_summary, dabg_summary_file, row.names = FALSE, quote = TRUE)

fit <- limma::lmFit(expression_matrix, design)
fit_contrasts <- limma::contrasts.fit(fit, contrast_matrix)
fit_ebayes <- limma::eBayes(fit_contrasts)

light_filter_keep <- dabg_detection_summary$dabg_light_filter_keep
filtered_expression_matrix <- expression_matrix[light_filter_keep, , drop = FALSE]
filtered_dabg_summary <- dabg_detection_summary[light_filter_keep, , drop = FALSE]

fit_filtered <- limma::lmFit(filtered_expression_matrix, design)
fit_filtered_contrasts <- limma::contrasts.fit(fit_filtered, contrast_matrix)
fit_filtered_ebayes <- limma::eBayes(fit_filtered_contrasts)

all_probe_results <- list()
filtered_results <- list()

for (contrast_name in colnames(contrast_matrix)) {
  all_output_file <- file.path(
    output_dir,
    paste0(
      "first_pass_limma_all_probes_",
      contrast_name,
      "_mirna_150001.csv"
    )
  )
  filtered_output_file <- file.path(
    output_dir,
    paste0(
      "first_pass_limma_dabg20_",
      contrast_name,
      "_mirna_150001.csv"
    )
  )

  all_probe_results[[contrast_name]] <- write_limma_results(
    fit_ebayes,
    contrast_name,
    all_output_file,
    dabg_detection_summary
  )
  filtered_results[[contrast_name]] <- write_limma_results(
    fit_filtered_ebayes,
    contrast_name,
    filtered_output_file,
    filtered_dabg_summary
  )
}

count_significant <- function(result, adj_p_cutoff = 0.05) {
  sum(result$adj.P.Val < adj_p_cutoff, na.rm = TRUE)
}

summary_lines <- c(
  "miRNA 150001 first-pass limma differential expression summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Scope:",
  "- Treatment-only limma model with no covariates.",
  "- Included treatments: Sham, MCAO1hr, MCAO3hr.",
  "- Excluded controls, N/A treatment rows, unmatched rows, and MCAO24hr because n=1.",
  "- QC-flagged samples were retained.",
  "- No pathway analysis or interpretation was run.",
  "",
  paste("Expression matrix:", expression_file),
  paste("DABG p-value matrix:", dabg_file),
  paste("Sample sheet:", sample_sheet_file),
  "",
  paste("Samples included:", ncol(expression_matrix)),
  paste("All probesets tested:", nrow(expression_matrix)),
  paste("DABG-light-filtered probesets tested:", nrow(filtered_expression_matrix)),
  paste(
    "DABG light filter:",
    "DABG p <",
    dabg_detection_p,
    "in at least",
    paste0(100 * dabg_detection_fraction, "%"),
    "of included samples"
  ),
  "",
  "Treatment counts:",
  capture.output(print(table(analysis_samples$treatment, useNA = "ifany"))),
  "",
  "Contrasts:",
  paste(names(contrast_definitions), contrast_definitions, sep = " = "),
  "",
  "Significant probesets at BH adjusted p < 0.05:",
  unlist(lapply(
    names(all_probe_results),
    function(contrast_name) {
      paste(
        contrast_name,
        "all_probes:",
        count_significant(all_probe_results[[contrast_name]])
      )
    }
  )),
  unlist(lapply(
    names(filtered_results),
    function(contrast_name) {
      paste(
        contrast_name,
        "dabg20:",
        count_significant(filtered_results[[contrast_name]])
      )
    }
  )),
  "",
  paste("Sample table:", sample_table_file),
  paste("Design matrix:", design_file),
  paste("Contrast matrix:", contrast_file),
  paste("DABG detection summary:", dabg_summary_file)
)

writeLines(summary_lines, summary_file)

message("Wrote first-pass limma sample table: ", sample_table_file)
message("Wrote first-pass limma design matrix: ", design_file)
message("Wrote first-pass limma contrast matrix: ", contrast_file)
message("Wrote first-pass DABG detection summary: ", dabg_summary_file)
message("Wrote first-pass limma summary: ", summary_file)
message("Wrote first-pass limma result tables to: ", output_dir)
