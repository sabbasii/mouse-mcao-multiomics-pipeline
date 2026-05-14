#!/usr/bin/env Rscript

# Run exploratory differential-expression analysis for the analysis-ready
# Clariom S Mouse HT mRNA dataset.
#
# Model argument:
#   treatment_only      ~ 0 + treatment
#   treatment_sex       ~ 0 + treatment + sex
#   treatment_age       ~ 0 + treatment + age_group
#   treatment_sex_age   ~ 0 + treatment + sex + age_group
#
# Inputs:
#   results/mrna/analysis_ready/
#     expression_matrix_unique_gene_mapped_mrna.csv
#       RMA expression for uniquely gene-mapped transcript clusters.
#     transcript_cluster_annotation_unique_gene_mapped_mrna.csv
#       One biological annotation row per expression transcript cluster.
#     analysis_samples_mrna.csv
#       Validated mRNA sample metadata used to confirm manifest assignments.
#
#   results/multiomics/sample_manifest/
#     paired_manifest_mirna_mrna.csv
#       The 43 paired animals eligible for analysis, including the exact mRNA
#       filename and standardized treatment for each animal.
#
# Contrasts:
#   MCAO1hr_vs_Sham      = MCAO1hr - Sham
#   MCAO3hr_vs_Sham      = MCAO3hr - Sham
#   MCAO3hr_vs_MCAO1hr   = MCAO3hr - MCAO1hr
#
# A positive logFC means higher expression in the first treatment named in a
# contrast; a negative logFC means higher expression in the second treatment.
#
# Outputs:
#   results/mrna/differential_expression/<model>_limma/
#     analysis_samples_mrna.csv
#       The 43 included samples in paired-manifest order.
#     design_matrix_mrna.csv
#       Design matrix used by the selected limma model.
#     contrast_matrix_mrna.csv
#       The three requested contrast definitions.
#     transcript_cluster_<contrast>_mrna.csv
#       Complete limma results for all retained transcript clusters, joined to
#       Entrez ID, gene symbol, gene name, gene type, and Ensembl ID.
#     limma_summary_mrna.txt
#       Model, sample counts, feature count, and raw/FDR result counts for each
#       contrast.
#
# Usage:
#   From the repository root:
#   Rscript scripts/mrna/05_differential_expression_mrna.R
#   Rscript scripts/mrna/05_differential_expression_mrna.R treatment_sex
#   Rscript scripts/mrna/05_differential_expression_mrna.R treatment_age
#   Rscript scripts/mrna/05_differential_expression_mrna.R treatment_sex_age

args <- commandArgs(trailingOnly = TRUE)
valid_models <- c(
  "treatment_only", "treatment_sex", "treatment_age",
  "treatment_sex_age"
)
if (length(args) > 1L) {
  stop(
    "Supply zero or one model argument: ",
    paste(valid_models, collapse = ", "),
    call. = FALSE
  )
}
analysis_model <- if (length(args) == 1L) {
  args[[1]]
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
  "results", "mrna", "analysis_ready",
  "expression_matrix_unique_gene_mapped_mrna.csv"
)
annotation_file <- file.path(
  "results", "mrna", "analysis_ready",
  "transcript_cluster_annotation_unique_gene_mapped_mrna.csv"
)
sample_file <- file.path(
  "results", "mrna", "analysis_ready", "analysis_samples_mrna.csv"
)
paired_manifest_file <- file.path(
  "results", "multiomics", "sample_manifest",
  "paired_manifest_mirna_mrna.csv"
)
output_dir <- file.path(
  "results", "mrna", "differential_expression",
  paste0(analysis_model, "_limma")
)

required_inputs <- c(
  expression_file,
  annotation_file,
  sample_file,
  paired_manifest_file
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Missing required package: limma", call. = FALSE)
}

expression_table <- read.csv(
  expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
annotation <- read.csv(
  annotation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
mrna_samples <- read.csv(
  sample_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)
paired_manifest <- read.csv(
  paired_manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)

required_expression_id <- "transcript_cluster_id"
if (!required_expression_id %in% names(expression_table)) {
  stop(
    "Expression matrix lacks transcript_cluster_id.",
    call. = FALSE
  )
}
if (!required_expression_id %in% names(annotation)) {
  stop(
    "Annotation table lacks transcript_cluster_id.",
    call. = FALSE
  )
}

required_manifest_columns <- c(
  "animal_id", "treatment", "sex", "age_group", "mrna_file_name"
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

required_sample_columns <- c(
  "file_name", "animal_id", "sex", "age_group", "treatment",
  "is_control", "match_status"
)
missing_sample_columns <- setdiff(
  required_sample_columns,
  names(mrna_samples)
)
if (length(missing_sample_columns) > 0L) {
  stop(
    "mRNA sample metadata lacks required column(s): ",
    paste(missing_sample_columns, collapse = ", "),
    call. = FALSE
  )
}

included_treatments <- c("Sham", "MCAO1hr", "MCAO3hr")
if (
  nrow(paired_manifest) != 43L ||
    anyDuplicated(paired_manifest$animal_id) ||
    anyDuplicated(paired_manifest$mrna_file_name) ||
    !all(paired_manifest$treatment %in% included_treatments)
) {
  stop(
    "Paired manifest must contain 43 unique Sham, MCAO1hr, and MCAO3hr mice.",
    call. = FALSE
  )
}

sample_index <- match(
  paired_manifest$mrna_file_name,
  mrna_samples$file_name
)
if (anyNA(sample_index)) {
  stop(
    "At least one paired mRNA filename is absent from the sample metadata.",
    call. = FALSE
  )
}
samples <- mrna_samples[sample_index, , drop = FALSE]

if (
  any(samples$is_control) ||
    any(samples$match_status != "matched") ||
    !identical(samples$animal_id, paired_manifest$animal_id) ||
    !identical(samples$treatment, paired_manifest$treatment) ||
    !identical(samples$sex, paired_manifest$sex) ||
    !identical(samples$age_group, paired_manifest$age_group)
) {
  stop(
    "Paired manifest and mRNA sample metadata do not agree.",
    call. = FALSE
  )
}
if (!all(samples$file_name %in% names(expression_table))) {
  stop(
    "At least one paired mRNA filename is absent from the expression matrix.",
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

if (
  any(!is.finite(expression)) ||
    anyDuplicated(rownames(expression)) ||
    anyDuplicated(annotation$transcript_cluster_id) ||
    !identical(
      rownames(expression),
      as.character(annotation$transcript_cluster_id)
    )
) {
  stop(
    "Expression or annotation validation failed.",
    call. = FALSE
  )
}

samples$treatment <- factor(
  samples$treatment,
  levels = included_treatments
)
samples$sex <- factor(samples$sex, levels = c("Female", "Male"))
samples$age_group <- factor(
  samples$age_group,
  levels = c("Old", "Young")
)
group_counts <- table(samples$treatment)
expected_group_counts <- c(
  Sham = 16L,
  MCAO1hr = 16L,
  MCAO3hr = 11L
)
if (!identical(
  unname(as.integer(group_counts)),
  unname(as.integer(expected_group_counts))
)) {
  stop(
    "Unexpected treatment counts in paired mRNA samples.",
    call. = FALSE
  )
}
if (analysis_model %in% c("treatment_sex", "treatment_sex_age")) {
  if (any(table(samples$sex) == 0L)) {
    stop("At least one required sex group has no samples.", call. = FALSE)
  }
}
if (analysis_model %in% c("treatment_age", "treatment_sex_age")) {
  if (any(table(samples$age_group) == 0L)) {
    stop("At least one required age group has no samples.", call. = FALSE)
  }
}

if (analysis_model == "treatment_only") {
  design <- stats::model.matrix(~ 0 + treatment, data = samples)
  model_description <- "treatment only (~0 + treatment)"
} else if (analysis_model == "treatment_sex") {
  design <- stats::model.matrix(
    ~ 0 + treatment + sex,
    data = samples
  )
  model_description <- "treatment + sex (~0 + treatment + sex)"
} else if (analysis_model == "treatment_age") {
  design <- stats::model.matrix(
    ~ 0 + treatment + age_group,
    data = samples
  )
  model_description <- paste(
    "treatment + age",
    "(~0 + treatment + age_group)"
  )
} else {
  design <- stats::model.matrix(
    ~ 0 + treatment + sex + age_group,
    data = samples
  )
  model_description <- paste(
    "treatment + sex + age",
    "(~0 + treatment + sex + age_group)"
  )
}
colnames(design) <- sub("^treatment", "", colnames(design))
if (qr(design)$rank != ncol(design)) {
  stop("The selected design matrix is not full rank.", call. = FALSE)
}

contrast_definitions <- c(
  MCAO1hr_vs_Sham = "MCAO1hr - Sham",
  MCAO3hr_vs_Sham = "MCAO3hr - Sham",
  MCAO3hr_vs_MCAO1hr = "MCAO3hr - MCAO1hr"
)
contrasts <- limma::makeContrasts(
  contrasts = contrast_definitions,
  levels = design
)

fit <- limma::lmFit(expression, design)
fit <- limma::contrasts.fit(fit, contrasts)
fit <- limma::eBayes(fit, robust = TRUE)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
write.csv(
  samples,
  file.path(output_dir, "analysis_samples_mrna.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
write.csv(
  design,
  file.path(output_dir, "design_matrix_mrna.csv"),
  row.names = TRUE,
  quote = TRUE
)
write.csv(
  contrasts,
  file.path(output_dir, "contrast_matrix_mrna.csv"),
  row.names = TRUE,
  quote = TRUE
)

result_counts <- character()
for (contrast_name in colnames(contrasts)) {
  result <- limma::topTable(
    fit,
    coef = contrast_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  result <- data.frame(
    transcript_cluster_id = rownames(result),
    result,
    row.names = NULL,
    check.names = FALSE
  )
  annotation_index <- match(
    result$transcript_cluster_id,
    annotation$transcript_cluster_id
  )
  if (anyNA(annotation_index)) {
    stop(
      "A limma result lacks its validated annotation row.",
      call. = FALSE
    )
  }
  result <- cbind(
    result,
    annotation[
      annotation_index,
      setdiff(names(annotation), "transcript_cluster_id"),
      drop = FALSE
    ]
  )
  write.csv(
    result,
    file.path(
      output_dir,
      paste0(
        "transcript_cluster_",
        contrast_name,
        "_mrna.csv"
      )
    ),
    row.names = FALSE,
    quote = TRUE,
    na = ""
  )
  result_counts <- c(
    result_counts,
    "",
    paste("Contrast:", contrast_name),
    paste("  Raw P < 0.05:  ", sum(result$P.Value < 0.05)),
    paste("  BH FDR < 0.05: ", sum(result$adj.P.Val < 0.05)),
    paste("  BH FDR < 0.10: ", sum(result$adj.P.Val < 0.10)),
    paste("  BH FDR < 0.20: ", sum(result$adj.P.Val < 0.20))
  )
}

writeLines(
  c(
    "Dataset: analysis-ready Clariom S Mouse HT mRNA",
    paste("Model:", model_description),
    "Interpretation status: exploratory, not final biology",
    paste("Paired manifest:", paired_manifest_file),
    paste("Samples:", ncol(expression)),
    paste("Transcript clusters:", nrow(expression)),
    paste(
      names(group_counts),
      as.integer(group_counts),
      collapse = "; "
    ),
    "Contrasts:",
    "  MCAO1hr_vs_Sham = MCAO1hr - Sham",
    "  MCAO3hr_vs_Sham = MCAO3hr - Sham",
    "  MCAO3hr_vs_MCAO1hr = MCAO3hr - MCAO1hr",
    result_counts
  ),
  file.path(output_dir, "limma_summary_mrna.txt")
)

message("Completed mRNA limma model: ", analysis_model)
message("Output: ", output_dir)
