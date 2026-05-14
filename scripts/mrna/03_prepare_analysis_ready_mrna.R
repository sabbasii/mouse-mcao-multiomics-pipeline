#!/usr/bin/env Rscript

# Prepare the primary analysis-ready Clariom S/mRNA expression dataset.
#
# Feature rule:
#   Keep non-control transcript clusters that map uniquely to one gene.
#   Exclude transcript clusters that are controls, absent from the biological
#   annotation database, lack a gene assignment, or map to multiple genes.
#
# Sample rule:
#   Keep biological samples with match_status == "matched".
#   Exclude control arrays.
#
# Prerequisites:
#   Rscript scripts/mrna/00_build_sample_sheet_mrna.R
#   Rscript scripts/mrna/01_array_qc_normalize_mrna.R
#   Rscript scripts/mrna/02_annotate_transcript_clusters_mrna.R
#
# Outputs:
#   results/mrna/analysis_ready/
#     expression_matrix_unique_gene_mapped_mrna.csv
#       Expression values for retained transcript clusters and biological
#       samples.
#     transcript_cluster_annotation_unique_gene_mapped_mrna.csv
#       One gene annotation row for each retained transcript cluster.
#     analysis_samples_mrna.csv
#       Sample metadata for the biological arrays retained in the matrix.
#     transcript_cluster_filter_report_mrna.csv
#       One row per input transcript cluster with its inclusion/exclusion reason.
#     analysis_ready_summary_mrna.txt
#       Filtering counts, dimensions, input paths, and validation results.
#
# Run from the repository root:
#   Rscript scripts/mrna/03_prepare_analysis_ready_mrna.R

expression_file <- file.path(
  "results", "mrna", "normalized",
  "rma_expression_mrna_150002.csv"
)
annotation_file <- file.path(
  "results", "mrna", "annotation",
  "transcript_cluster_gene_annotation_mrna_150002.csv"
)
sample_sheet_file <- file.path(
  "results", "mrna", "sample_sheet",
  "sample_sheet_mrna_150002.csv"
)
output_dir <- file.path("results", "mrna", "analysis_ready")

analysis_expression_file <- file.path(
  output_dir,
  "expression_matrix_unique_gene_mapped_mrna.csv"
)
analysis_annotation_file <- file.path(
  output_dir,
  "transcript_cluster_annotation_unique_gene_mapped_mrna.csv"
)
analysis_samples_file <- file.path(
  output_dir,
  "analysis_samples_mrna.csv"
)
filter_report_file <- file.path(
  output_dir,
  "transcript_cluster_filter_report_mrna.csv"
)
summary_file <- file.path(
  output_dir,
  "analysis_ready_summary_mrna.txt"
)

required_inputs <- c(
  expression_file,
  annotation_file,
  sample_sheet_file
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

expression <- read.csv(
  expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
annotation <- read.csv(
  annotation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
sample_sheet <- read.csv(
  sample_sheet_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"transcript_cluster_id" %in% names(expression)) {
  stop("Expression matrix lacks transcript_cluster_id.", call. = FALSE)
}
if (!all(c(
  "transcript_cluster_id",
  "mapping_status",
  "gene_mapping_count",
  "ENTREZID",
  "SYMBOL",
  "GENENAME",
  "GENETYPE",
  "ENSEMBL"
) %in% names(annotation))) {
  stop("Annotation table lacks required fields.", call. = FALSE)
}
if (!all(c(
  "file_name",
  "animal_id",
  "is_control",
  "match_status"
) %in% names(sample_sheet))) {
  stop("Sample sheet lacks required fields.", call. = FALSE)
}

expression_ids <- as.character(expression$transcript_cluster_id)
if (anyDuplicated(expression_ids)) {
  stop("Expression matrix contains duplicate transcript-cluster IDs.", call. = FALSE)
}
if (!setequal(expression_ids, unique(annotation$transcript_cluster_id))) {
  stop(
    "Expression and annotation transcript-cluster ID sets differ.",
    call. = FALSE
  )
}

biological_samples <- sample_sheet[
  sample_sheet$match_status == "matched" &
    !sample_sheet$is_control,
  ,
  drop = FALSE
]
sample_columns <- biological_samples$file_name
if (anyDuplicated(sample_columns)) {
  stop("Biological sample filenames are duplicated.", call. = FALSE)
}
if (!all(sample_columns %in% names(expression))) {
  stop(
    "Not every biological sample is present in the expression matrix.",
    call. = FALSE
  )
}

feature_status <- unique(annotation[, c(
  "transcript_cluster_id",
  "mapping_status",
  "gene_mapping_count"
)])
status_count_per_id <- table(feature_status$transcript_cluster_id)
if (any(status_count_per_id != 1L)) {
  stop(
    "A transcript cluster has inconsistent mapping status or count.",
    call. = FALSE
  )
}
feature_status <- feature_status[
  match(expression_ids, feature_status$transcript_cluster_id),
  ,
  drop = FALSE
]

is_control_feature <- grepl(
  "^AFFX",
  feature_status$transcript_cluster_id,
  ignore.case = TRUE
)
feature_status$filter_status <- ifelse(
  is_control_feature,
  "excluded_control_feature",
  ifelse(
    feature_status$mapping_status == "not_in_annotation_database",
    "excluded_not_in_annotation_database",
    ifelse(
      feature_status$mapping_status == "database_without_gene",
      "excluded_without_gene_assignment",
      ifelse(
        feature_status$gene_mapping_count > 1L,
        "excluded_one_to_many_gene_mapping",
        "retained_unique_gene_mapping"
      )
    )
  )
)

retained_ids <- feature_status$transcript_cluster_id[
  feature_status$filter_status == "retained_unique_gene_mapping"
]
if (length(retained_ids) == 0L) {
  stop("Feature filtering retained no transcript clusters.", call. = FALSE)
}

retained_annotation_long <- annotation[
  annotation$transcript_cluster_id %in% retained_ids,
  ,
  drop = FALSE
]
core_annotation <- unique(retained_annotation_long[, c(
  "transcript_cluster_id",
  "ENTREZID",
  "SYMBOL",
  "GENENAME",
  "GENETYPE"
)])
if (anyDuplicated(core_annotation$transcript_cluster_id)) {
  stop(
    "A supposedly unique gene mapping has conflicting core gene annotations.",
    call. = FALSE
  )
}

collapse_identifiers <- function(values) {
  values <- sort(unique(values[!is.na(values) & nzchar(values)]))
  if (length(values) == 0L) {
    return(NA_character_)
  }
  paste(values, collapse = " /// ")
}
ensembl_by_transcript_cluster <- tapply(
  retained_annotation_long$ENSEMBL,
  retained_annotation_long$transcript_cluster_id,
  collapse_identifiers
)
core_annotation$ENSEMBL <- unname(
  ensembl_by_transcript_cluster[core_annotation$transcript_cluster_id]
)

analysis_annotation <- core_annotation[
  match(retained_ids, core_annotation$transcript_cluster_id),
  ,
  drop = FALSE
]
if (
  anyNA(analysis_annotation$transcript_cluster_id) ||
    anyDuplicated(analysis_annotation$transcript_cluster_id)
) {
  stop("Analysis annotation is incomplete or duplicated.", call. = FALSE)
}

analysis_expression <- expression[
  match(retained_ids, expression_ids),
  c("transcript_cluster_id", sample_columns),
  drop = FALSE
]
expression_values <- as.matrix(
  analysis_expression[, sample_columns, drop = FALSE]
)
storage.mode(expression_values) <- "numeric"
if (any(!is.finite(expression_values))) {
  stop("Analysis-ready expression contains non-finite values.", call. = FALSE)
}
if (!identical(
  analysis_expression$transcript_cluster_id,
  analysis_annotation$transcript_cluster_id
)) {
  stop(
    "Analysis-ready expression and annotation orders differ.",
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  analysis_expression,
  analysis_expression_file,
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  analysis_annotation,
  analysis_annotation_file,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
write.csv(
  biological_samples,
  analysis_samples_file,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
write.csv(
  feature_status,
  filter_report_file,
  row.names = FALSE,
  quote = TRUE
)

filter_counts <- table(feature_status$filter_status)
writeLines(
  c(
    "Dataset: Clariom S Mouse HT mRNA 150002",
    paste("Expression input:", expression_file),
    paste("Annotation input:", annotation_file),
    paste("Sample-sheet input:", sample_sheet_file),
    "",
    "Primary feature rule:",
    "Retain non-control transcript clusters with exactly one gene mapping.",
    "Exclude controls, database-unmatched IDs, gene-less IDs, and one-to-many mappings.",
    "",
    paste("Input transcript clusters:", length(expression_ids)),
    paste("Retained transcript clusters:", length(retained_ids)),
    paste("Input arrays:", ncol(expression) - 1L),
    paste("Retained biological samples:", length(sample_columns)),
    paste(
      "Excluded control arrays:",
      sum(sample_sheet$is_control)
    ),
    "",
    "Feature filter counts:",
    paste(
      paste(names(filter_counts), as.integer(filter_counts), sep = ": "),
      collapse = "\n"
    ),
    "",
    "Validation:",
    "Expression and annotation transcript-cluster order is identical: TRUE",
    "Duplicate retained transcript-cluster IDs: 0",
    "Non-finite analysis-ready expression values: 0"
  ),
  summary_file
)

message("Wrote analysis-ready expression: ", analysis_expression_file)
message("Wrote analysis-ready annotation: ", analysis_annotation_file)
message("Wrote analysis samples: ", analysis_samples_file)
message("Wrote feature-filter report: ", filter_report_file)
message("Wrote summary: ", summary_file)
