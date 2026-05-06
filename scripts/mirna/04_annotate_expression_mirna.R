#!/usr/bin/env Rscript

# Annotate and filter a prepared miRNA expression dataset before differential
# expression analysis.
#
# This script:
#   1. matches ProbeSetName values in the expression matrix to Thermo's
#      miRNA 4.1 annotation using the "Probe Set ID" field;
#   2. retains probes annotated to Mus musculus;
#   3. retains the mature-miRNA feature class ("Sequence Type" == "miRNA");
#   4. writes the matched annotation, the filtered mouse mature-miRNA
#      annotation, and a DE-ready expression matrix.
#
# Required argument:
#   dataset  One of: vendor_chp, complete_rma_excluding_two, rma_normalized_mirna
#
# Outputs:
#   results/mirna/expression/<dataset>/annotation/
#     matched_probeset_annotation_mirna.csv
#     mouse_mature_mirna_annotation_mirna.csv
#     mouse_mature_mirna_expression_mirna.csv
#     annotation_summary_mirna.txt
#
# Usage:
#   Rscript scripts/mirna/04_annotate_expression_mirna.R vendor_chp
#   Rscript scripts/mirna/04_annotate_expression_mirna.R complete_rma_excluding_two
#   Rscript scripts/mirna/04_annotate_expression_mirna.R rma_normalized_mirna

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

expression_dir <- file.path("results", "mirna", "expression", dataset)
expression_file <- file.path(
  expression_dir, "expression_matrix_mirna.csv"
)
annotation_file <- file.path(
  "resources", "affymetrix_library_files", "mirna_4_1",
  "miRNA-4_1-st-v1.annotations.20160922.csv"
)
output_dir <- file.path(expression_dir, "annotation")
matched_annotation_file <- file.path(
  output_dir, "matched_probeset_annotation_mirna.csv"
)
filtered_annotation_file <- file.path(
  output_dir, "mouse_mature_mirna_annotation_mirna.csv"
)
filtered_expression_file <- file.path(
  output_dir, "mouse_mature_mirna_expression_mirna.csv"
)
summary_file <- file.path(
  output_dir, "annotation_summary_mirna.txt"
)

for (path in c(expression_file, annotation_file)) {
  if (!file.exists(path)) {
    stop("Required input not found: ", path, call. = FALSE)
  }
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
  comment.char = "#"
)

required_expression_columns <- "ProbeSetName"
required_annotation_columns <- c(
  "Probe Set ID",
  "Probe Set Name",
  "Accession",
  "Transcript ID(Array Design)",
  "Sequence Type",
  "Species Scientific Name",
  "Sequence",
  "Annotation Date",
  "Sequence Source"
)
if (!all(required_expression_columns %in% names(expression))) {
  stop("Expression matrix lacks ProbeSetName.", call. = FALSE)
}
if (!all(required_annotation_columns %in% names(annotation))) {
  stop("Thermo annotation file lacks required columns.", call. = FALSE)
}

expression$ProbeSetName <- as.character(expression$ProbeSetName)
annotation[["Probe Set ID"]] <- as.character(annotation[["Probe Set ID"]])

if (anyDuplicated(expression$ProbeSetName)) {
  stop("Expression matrix contains duplicated ProbeSetName values.", call. = FALSE)
}
if (anyDuplicated(annotation[["Probe Set ID"]])) {
  stop("Thermo annotation contains duplicated Probe Set ID values.", call. = FALSE)
}

annotation_index <- match(
  expression$ProbeSetName,
  annotation[["Probe Set ID"]]
)
matched_count <- sum(!is.na(annotation_index))
unmatched_count <- sum(is.na(annotation_index))
if (matched_count == 0L) {
  stop("No expression probesets matched the Thermo annotation.", call. = FALSE)
}

selected_annotation <- annotation[
  annotation_index,
  required_annotation_columns,
  drop = FALSE
]
names(selected_annotation)[names(selected_annotation) == "Probe Set ID"] <-
  "ProbeSetName"

matched_annotation <- selected_annotation[!is.na(annotation_index), , drop = FALSE]
mouse_mature_keep <- (
  !is.na(annotation_index) &
    selected_annotation[["Species Scientific Name"]] == "Mus musculus" &
    selected_annotation[["Sequence Type"]] == "miRNA"
)

if (!any(mouse_mature_keep)) {
  stop("No Mus musculus mature-miRNA probesets were identified.", call. = FALSE)
}

mouse_mature_annotation <- selected_annotation[
  mouse_mature_keep, , drop = FALSE
]
mouse_mature_expression <- expression[
  mouse_mature_keep, , drop = FALSE
]

if (!identical(
  mouse_mature_expression$ProbeSetName,
  mouse_mature_annotation$ProbeSetName
)) {
  stop("Filtered expression and annotation orders do not match.", call. = FALSE)
}

sample_columns <- setdiff(
  names(mouse_mature_expression),
  c("ProbeSetName", "feature_name", "ID")
)
expression_values <- as.matrix(
  mouse_mature_expression[, sample_columns, drop = FALSE]
)
storage.mode(expression_values) <- "numeric"
if (any(!is.finite(expression_values))) {
  stop("Filtered expression matrix contains non-finite values.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  matched_annotation,
  matched_annotation_file,
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  mouse_mature_annotation,
  filtered_annotation_file,
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  mouse_mature_expression,
  filtered_expression_file,
  row.names = FALSE,
  quote = TRUE
)

writeLines(
  c(
    paste("Dataset:", dataset),
    paste("Expression input:", expression_file),
    paste("Thermo annotation input:", annotation_file),
    "Annotation release: September 22, 2016",
    "miRBase version: v20",
    paste("Expression probesets:", nrow(expression)),
    paste("Matched to annotation:", matched_count),
    paste("Unmatched:", unmatched_count),
    "DE-ready filter: Species Scientific Name == Mus musculus",
    "DE-ready filter: Sequence Type == miRNA",
    paste("Mouse mature-miRNA probesets retained:", nrow(mouse_mature_expression)),
    paste("Samples retained:", length(sample_columns))
  ),
  summary_file
)

message("Annotated dataset: ", dataset)
message("DE-ready expression matrix: ", filtered_expression_file)
