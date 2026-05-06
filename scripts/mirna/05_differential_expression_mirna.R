#!/usr/bin/env Rscript

# Run exploratory differential-expression analysis for one annotated miRNA
# expression dataset.
#
# This script:
#   1. reads the DE-ready Mus musculus mature-miRNA expression matrix and its
#      matching Thermo annotation created by script 04;
#   2. reads the paired miRNA-mRNA manifest and selects its animals directly
#      from the existing expression matrix without writing a second matrix;
#   3. matches manifest filenames to the miRNA sample sheet and preserves the
#      manifest animal order;
#   4. fits a treatment-only model or an additive model adjusting for sex or
#      age;
#   5. tests MCAO1hr versus Sham, MCAO3hr versus Sham, and MCAO3hr
#      versus MCAO1hr;
#   6. analyzes all retained mouse mature-miRNA probesets and, when vendor
#      DABG values are available, a DABG20-filtered subset;
#   7. adds miRNA names, accessions, sequence information, and detection
#      summaries to the result tables.
#
# Arguments:
#   dataset  One of: vendor_chp, complete_rma_excluding_two, rma_normalized_mirna
#   model    Optional; one of: treatment_only, treatment_sex, treatment_age,
#            treatment_sex_age
#            Default: treatment_only
#
# Required prior step:
#   Rscript scripts/mirna/04_annotate_expression_mirna.R <dataset>
#
# Outputs:
#   results/mirna/differential_expression/<dataset>/<model>_limma/
#     analysis_samples_mirna.csv
#     design_matrix_mirna.csv
#     contrast_matrix_mirna.csv
#     mouse_mature_mirna_<contrast>_mirna.csv
#     dabg20_<contrast>_mirna.csv
#     dabg_detection_summary_mirna.csv
#     limma_summary_mirna.txt
#
# Usage:
  # Rscript scripts/mirna/05_differential_expression_mirna.R vendor_chp
  # Rscript scripts/mirna/05_differential_expression_mirna.R complete_rma_excluding_two
  # Rscript scripts/mirna/05_differential_expression_mirna.R rma_normalized_mirna
  # Rscript scripts/mirna/05_differential_expression_mirna.R rma_normalized_mirna treatment_sex
  # Rscript scripts/mirna/05_differential_expression_mirna.R rma_normalized_mirna treatment_age
  # Rscript scripts/mirna/05_differential_expression_mirna.R rma_normalized_mirna treatment_sex_age

args <- commandArgs(trailingOnly = TRUE)
valid_datasets <- c(
  "vendor_chp", "complete_rma_excluding_two", "rma_normalized_mirna"
)
valid_models <- c(
  "treatment_only", "treatment_sex", "treatment_age",
  "treatment_sex_age"
)
if (
  length(args) < 1L ||
    length(args) > 2L ||
    !args[[1]] %in% valid_datasets
) {
  stop(
    "Usage: Rscript scripts/mirna/05_differential_expression_mirna.R ",
    "<dataset> [treatment_only|treatment_sex|treatment_age|",
    "treatment_sex_age]. Valid datasets: ",
    paste(valid_datasets, collapse = ", "),
    call. = FALSE
  )
}
dataset <- args[[1]]
analysis_model <- if (length(args) == 2L) {
  args[[2]]
} else {
  "treatment_only"
}
if (!analysis_model %in% valid_models) {
  stop(
    "Unknown model: ",
    analysis_model,
    ". Choose: ",
    paste(valid_models, collapse = ", "),
    call. = FALSE
  )
}

expression_file <- file.path(
  "results", "mirna", "expression", dataset,
  "annotation", "mouse_mature_mirna_expression_mirna.csv"
)
annotation_file <- file.path(
  "results", "mirna", "expression", dataset,
  "annotation", "mouse_mature_mirna_annotation_mirna.csv"
)
dabg_file <- file.path(
  "results", "mirna", "expression", "vendor_chp",
  "dabg_pvalues_mirna.csv"
)
sample_sheet_file <- file.path(
  "results", "mirna", "sample_sheet", "sample_sheet_mirna.csv"
)
paired_manifest_file <- file.path(
  "results", "multiomics", "sample_manifest",
  "paired_manifest_mirna_mrna.csv"
)
output_dir <- file.path(
  "results", "mirna", "differential_expression", dataset,
  paste0(analysis_model, "_limma")
)
included_treatments <- c("Sham", "MCAO1hr", "MCAO3hr")
contrast_definitions <- c(
  MCAO1hr_vs_Sham = "MCAO1hr - Sham",
  MCAO3hr_vs_Sham = "MCAO3hr - Sham",
  MCAO3hr_vs_MCAO1hr = "MCAO3hr - MCAO1hr"
)

read_probe_matrix <- function(path) {
  table <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
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

write_results <- function(
  fit, coefficient, path, annotation_table, detection = NULL
) {
  result <- limma::topTable(
    fit, coef = coefficient, number = Inf,
    adjust.method = "BH", sort.by = "P"
  )
  result <- data.frame(
    ProbeSetName = rownames(result), result,
    row.names = NULL, check.names = FALSE
  )
  result <- cbind(
    result,
    annotation_table[
      match(result$ProbeSetName, annotation_table$ProbeSetName),
      setdiff(names(annotation_table), "ProbeSetName"),
      drop = FALSE
    ]
  )
  if (!is.null(detection)) {
    result <- cbind(
      result,
      detection[match(result$ProbeSetName, detection$ProbeSetName), -1, drop = FALSE]
    )
  }
  write.csv(result, path, row.names = FALSE, quote = TRUE)
  result
}

if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Missing required package: limma", call. = FALSE)
}
if (
  !file.exists(expression_file) ||
    !file.exists(annotation_file) ||
    !file.exists(sample_sheet_file) ||
    !file.exists(paired_manifest_file)
) {
  stop(
    "Required annotated expression matrix, annotation, sample sheet, or paired ",
    "manifest is missing. Run script 04 and the paired-manifest step first.",
    call. = FALSE
  )
}

expression <- read_probe_matrix(expression_file)
if (any(!is.finite(expression))) {
  stop("Expression matrix contains non-finite values.", call. = FALSE)
}
annotation <- read.csv(
  annotation_file, stringsAsFactors = FALSE, check.names = FALSE
)
annotation$ProbeSetName <- as.character(annotation$ProbeSetName)
if (
  anyDuplicated(annotation$ProbeSetName) ||
    !setequal(rownames(expression), annotation$ProbeSetName)
) {
  stop("DE-ready expression and annotation probesets do not agree.", call. = FALSE)
}
sample_sheet <- read.csv(
  sample_sheet_file, stringsAsFactors = FALSE, check.names = FALSE
)
sample_sheet$sample_name <- sub(
  "[.]CEL$", "", sample_sheet$file_name, ignore.case = TRUE
)
paired_manifest <- read.csv(
  paired_manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_manifest_columns <- c(
  "animal_id", "treatment", "sex", "age_group", "mirna_file_name"
)
missing_manifest_columns <- setdiff(
  required_manifest_columns,
  names(paired_manifest)
)
if (length(missing_manifest_columns) > 0L) {
  stop(
    "Paired manifest lacks required column(s): ",
    paste(missing_manifest_columns, collapse = ", "),
    call. = FALSE
  )
}
if (
  anyDuplicated(paired_manifest$animal_id) ||
    anyDuplicated(paired_manifest$mirna_file_name) ||
    !all(paired_manifest$treatment %in% included_treatments)
) {
  stop(
    "Paired manifest contains duplicated samples or unexpected treatments.",
    call. = FALSE
  )
}

paired_manifest$sample_name <- sub(
  "[.]CEL$", "", paired_manifest$mirna_file_name, ignore.case = TRUE
)
available_manifest <- paired_manifest[
  paired_manifest$sample_name %in% colnames(expression),
  ,
  drop = FALSE
]
missing_paired_samples <- paired_manifest[
  !paired_manifest$sample_name %in% colnames(expression),
  ,
  drop = FALSE
]
sample_index <- match(
  available_manifest$mirna_file_name,
  sample_sheet$file_name
)
if (anyNA(sample_index)) {
  stop(
    "At least one paired-manifest miRNA filename is absent from the sample sheet.",
    call. = FALSE
  )
}
samples <- sample_sheet[sample_index, , drop = FALSE]
if (
  any(samples$is_control) ||
    any(samples$match_status != "matched") ||
    !identical(samples$animal_id, available_manifest$animal_id) ||
    !identical(samples$treatment, available_manifest$treatment) ||
    !identical(samples$sex, available_manifest$sex) ||
    !identical(samples$age, available_manifest$age_group)
) {
  stop(
    "Paired manifest and miRNA sample metadata do not agree.",
    call. = FALSE
  )
}
expression <- expression[, samples$sample_name, drop = FALSE]
samples$treatment <- factor(samples$treatment, levels = included_treatments)
samples$sex <- factor(samples$sex, levels = c("Female", "Male"))
samples$age <- factor(samples$age, levels = c("Old", "Young"))

if (any(table(samples$treatment) == 0L)) {
  stop("At least one required treatment group has no samples.", call. = FALSE)
}
if (analysis_model == "treatment_sex" && any(table(samples$sex) == 0L)) {
  stop("At least one required sex group has no samples.", call. = FALSE)
}
if (analysis_model == "treatment_age" && any(table(samples$age) == 0L)) {
  stop("At least one required age group has no samples.", call. = FALSE)
}
if (
  analysis_model == "treatment_sex_age" &&
    (any(table(samples$sex) == 0L) || any(table(samples$age) == 0L))
) {
  stop(
    "Both sex and age groups are required for treatment_sex_age.",
    call. = FALSE
  )
}

if (analysis_model == "treatment_only") {
  design <- model.matrix(~ 0 + treatment, data = samples)
  model_description <- "treatment only (~0 + treatment)"
} else if (analysis_model == "treatment_sex") {
  design <- model.matrix(~ 0 + treatment + sex, data = samples)
  model_description <- "treatment + sex (~0 + treatment + sex)"
} else if (analysis_model == "treatment_age") {
  design <- model.matrix(~ 0 + treatment + age, data = samples)
  model_description <- "treatment + age (~0 + treatment + age)"
} else {
  design <- model.matrix(
    ~ 0 + treatment + sex + age,
    data = samples
  )
  model_description <- paste(
    "treatment + sex + age",
    "(~0 + treatment + sex + age)"
  )
}
colnames(design) <- sub("^treatment", "", colnames(design))
if (qr(design)$rank != ncol(design)) {
  stop("The selected design matrix is not full rank.", call. = FALSE)
}
contrasts <- limma::makeContrasts(
  contrasts = contrast_definitions,
  levels = design
)

fit_expression <- function(matrix_object) {
  fit <- limma::lmFit(matrix_object, design)
  fit <- limma::contrasts.fit(fit, contrasts)
  limma::eBayes(fit, robust = TRUE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  samples,
  file.path(output_dir, "analysis_samples_mirna.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  design,
  file.path(output_dir, "design_matrix_mirna.csv"),
  row.names = TRUE,
  quote = TRUE
)
write.csv(
  contrasts,
  file.path(output_dir, "contrast_matrix_mirna.csv"),
  row.names = TRUE,
  quote = TRUE
)

all_fit <- fit_expression(expression)
result_counts <- character()
for (contrast_name in colnames(contrasts)) {
  result <- write_results(
    all_fit, contrast_name,
    file.path(
      output_dir,
      paste0("mouse_mature_mirna_", contrast_name, "_mirna.csv")
    ),
    annotation
  )
  result_counts <- c(
    result_counts,
    "",
    paste("Result set: mouse_mature_mirna | Contrast:", contrast_name),
    paste("  Raw P < 0.05:  ", sum(result$P.Value < 0.05)),
    paste("  BH FDR < 0.05: ", sum(result$adj.P.Val < 0.05)),
    paste("  BH FDR < 0.10: ", sum(result$adj.P.Val < 0.10)),
    paste("  BH FDR < 0.20: ", sum(result$adj.P.Val < 0.20))
  )
}

# The vendor DABG values are detection measurements, so the same light
# detection filter can be applied to shared probes/samples for any of the
# three expression normalizations.
if (file.exists(dabg_file)) {
  dabg <- read_probe_matrix(dabg_file)
  shared_samples <- samples$sample_name
  if (all(shared_samples %in% colnames(dabg))) {
    dabg <- dabg[, shared_samples, drop = FALSE]
    detection <- data.frame(
      ProbeSetName = rownames(dabg),
      dabg_detected_p_lt_0_05 = rowSums(dabg < 0.05),
      dabg_detected_fraction_p_lt_0_05 = rowMeans(dabg < 0.05),
      stringsAsFactors = FALSE
    )
    detection$dabg_light_filter_keep <-
      detection$dabg_detected_fraction_p_lt_0_05 >= 0.20
    write.csv(
      detection,
      file.path(output_dir, "dabg_detection_summary_mirna.csv"),
      row.names = FALSE,
      quote = TRUE
    )
    kept_ids <- detection$ProbeSetName[detection$dabg_light_filter_keep]
    kept_ids <- intersect(rownames(expression), kept_ids)
    filtered_fit <- fit_expression(expression[kept_ids, , drop = FALSE])
    for (contrast_name in colnames(contrasts)) {
      result <- write_results(
        filtered_fit, contrast_name,
        file.path(
          output_dir,
          paste0("dabg20_", contrast_name, "_mirna.csv")
        ), annotation,
        detection
      )
      result_counts <- c(
        result_counts,
        "",
        paste("Result set: dabg20 | Contrast:", contrast_name),
        paste("  Raw P < 0.05:  ", sum(result$P.Value < 0.05)),
        paste("  BH FDR < 0.05: ", sum(result$adj.P.Val < 0.05)),
        paste("  BH FDR < 0.10: ", sum(result$adj.P.Val < 0.10)),
        paste("  BH FDR < 0.20: ", sum(result$adj.P.Val < 0.20))
      )
    }
  }
}

writeLines(
  c(
    paste("Dataset:", dataset),
    paste("Model:", model_description),
    "Interpretation status: exploratory, not final biology",
    paste("Paired manifest:", paired_manifest_file),
    paste("Paired manifest samples:", nrow(paired_manifest)),
    paste("Samples:", ncol(expression)),
    paste(
      "Paired samples unavailable in this expression dataset:",
      nrow(missing_paired_samples)
    ),
    paste(
      "Unavailable paired animal IDs:",
      if (nrow(missing_paired_samples) == 0L) {
        "none"
      } else {
        paste(missing_paired_samples$animal_id, collapse = ", ")
      }
    ),
    "Feature filter: Mus musculus, Sequence Type == miRNA",
    paste("Mouse mature-miRNA probesets:", nrow(expression)),
    paste(names(table(samples$treatment)), table(samples$treatment), collapse = "; "),
    result_counts
  ),
  file.path(output_dir, "limma_summary_mirna.txt")
)

message("Completed exploratory limma for: ", dataset)
message("Output: ", output_dir)
