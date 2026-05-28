#!/usr/bin/env Rscript

# Combine standardized miRTarBase and TargetScanMouse relationships.
#
# Inputs:
#   results/multiomics/miRNA_target_evidence/mirtarbase/
#   results/multiomics/miRNA_target_evidence/targetscanmouse/
#
# Output:
#   results/multiomics/miRNA_target_evidence/combined/
#     mirna_mouse_target_evidence.csv
#
# Usage:
#   Rscript scripts/multiomics/03_combine_mirna_target_evidence.R

input_files <- c(
  mirtarbase = file.path(
    "results", "multiomics", "miRNA_target_evidence", "mirtarbase",
    "mirtarbase_mouse_targets_standardized.csv"
  ),
  targetscanmouse = file.path(
    "results", "multiomics", "miRNA_target_evidence", "targetscanmouse",
    "targetscanmouse_targets_standardized.csv"
  )
)
output_dir <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "combined"
)
output_file <- file.path(output_dir, "mirna_mouse_target_evidence.csv")

missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing standardized target table(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

read_standardized <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

mirtarbase <- read_standardized(input_files[["mirtarbase"]])
targetscanmouse <- read_standardized(input_files[["targetscanmouse"]])

common_columns <- union(names(mirtarbase), names(targetscanmouse))
add_missing_columns <- function(table, columns) {
  missing <- setdiff(columns, names(table))
  for (column in missing) table[[column]] <- NA
  table[, columns, drop = FALSE]
}

combined <- rbind(
  add_missing_columns(mirtarbase, common_columns),
  add_missing_columns(targetscanmouse, common_columns)
)

combined$evidence_sources <- combined$evidence_source
pair_key <- paste(combined$mirna_id, combined$target_gene_symbol, sep = "||")
source_by_pair <- tapply(
  combined$evidence_source,
  pair_key,
  function(x) paste(sort(unique(x)), collapse = ";")
)
combined$evidence_sources <- unname(source_by_pair[pair_key])

combined <- combined[order(combined$mirna_id, combined$target_gene_symbol), , drop = FALSE]
combined <- combined[!duplicated(combined), , drop = FALSE]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(combined, output_file, row.names = FALSE, na = "")

cat("Wrote ", nrow(combined), " combined target-evidence rows: ",
    output_file, "\n", sep = "")
