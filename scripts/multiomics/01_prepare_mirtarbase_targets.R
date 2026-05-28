#!/usr/bin/env Rscript

# Prepare the mouse miRTarBase 10.0 interaction table.
#
# Input:
#   resources/mirna_target_databases/mirtarbase_10.0/mmu_MTI.csv
#
# Output:
#   results/multiomics/miRNA_target_evidence/mirtarbase/
#     mirtarbase_mouse_targets_standardized.csv
#
# The output keeps one row per distinct mouse miRNA--gene relationship and
# retains the experimental evidence fields for later evidence ranking.
#
# Usage:
#   Rscript scripts/multiomics/01_prepare_mirtarbase_targets.R

input_file <- file.path(
  "resources", "mirna_target_databases", "mirtarbase_10.0", "mmu_MTI.csv"
)
output_dir <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "mirtarbase"
)
output_file <- file.path(
  output_dir, "mirtarbase_mouse_targets_standardized.csv"
)

if (!file.exists(input_file)) {
  stop("Missing miRTarBase input: ", input_file, call. = FALSE)
}

targets <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)

required_columns <- c(
  "miRNA", "Species (miRNA)", "Target Gene", "Target Gene (Entrez ID)",
  "Species (Target Gene)", "Experiments", "Support Type",
  "References (PMID)"
)
missing_columns <- setdiff(required_columns, names(targets))
if (length(missing_columns) > 0L) {
  stop(
    "miRTarBase is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

normalize_id <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.0+$", "", x)
}

mirna_ids <- tolower(trimws(targets[["miRNA"]]))
keep_mouse <-
  tolower(trimws(targets[["Species (miRNA)"]])) == "mmu" &
    tolower(trimws(targets[["Species (Target Gene)"]])) == "mmu" &
    grepl("^mmu-", mirna_ids)

targets <- targets[keep_mouse, , drop = FALSE]
mirna_ids <- mirna_ids[keep_mouse]

standardized <- data.frame(
  mirna_id = mirna_ids,
  target_gene_symbol = trimws(targets[["Target Gene"]]),
  target_entrez_id = normalize_id(targets[["Target Gene (Entrez ID)"]]),
  target_transcript_id = NA_character_,
  evidence_source = "miRTarBase",
  evidence_type = trimws(targets[["Support Type"]]),
  evidence_detail = trimws(targets[["Experiments"]]),
  reference_id = normalize_id(targets[["References (PMID)"]]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

standardized <- unique(standardized)
standardized <- standardized[
  nzchar(standardized$mirna_id) &
    nzchar(standardized$target_gene_symbol),
  ,
  drop = FALSE
]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(standardized, output_file, row.names = FALSE, na = "")

cat("Wrote ", nrow(standardized), " standardized miRTarBase relationships: ",
    output_file, "\n", sep = "")
