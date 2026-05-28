#!/usr/bin/env Rscript

# Prepare the mouse TargetScanMouse 8.0 default prediction table.
#
# Inputs:
#   resources/mirna_target_databases/targetscanmouse_8.0/
#     miR_Family_Info.txt
#     Gene_info.txt
#     Summary_Counts.default_predictions.txt
#
# Output:
#   results/multiomics/miRNA_target_evidence/targetscanmouse/
#     targetscanmouse_targets_standardized.csv
#
# Usage:
#   Rscript scripts/multiomics/02_prepare_targetscanmouse_targets.R

input_dir <- file.path(
  "resources", "mirna_target_databases", "targetscanmouse_8.0"
)
output_dir <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "targetscanmouse"
)
output_file <- file.path(
  output_dir, "targetscanmouse_targets_standardized.csv"
)

input_files <- file.path(
  input_dir,
  c(
    "miR_Family_Info.txt",
    "Gene_info.txt",
    "Summary_Counts.default_predictions.txt"
  )
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing TargetScanMouse input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

read_targetscan <- function(path) {
  read.delim(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "NULL")
  )
}

mirna_family <- read_targetscan(input_files[1])
gene_info <- read_targetscan(input_files[2])
predictions <- read_targetscan(input_files[3])

required_family <- c("miR family", "Species ID", "MiRBase ID")
required_gene <- c("Transcript ID", "Gene ID", "Gene symbol", "Species ID")
required_prediction <- c(
  "Transcript ID", "Gene Symbol", "miRNA family", "Species ID",
  "Representative miRNA", "Total context++ score",
  "Cumulative weighted context++ score", "Aggregate PCT"
)
for (required in list(required_family, required_gene, required_prediction)) {
  table_name <- if (identical(required, required_family)) {
    "miR_Family_Info.txt"
  } else if (identical(required, required_gene)) {
    "Gene_info.txt"
  } else {
    "Summary_Counts.default_predictions.txt"
  }
  missing_columns <- setdiff(required, names(
    list(mirna_family, gene_info, predictions)[[
      match(table_name, c("miR_Family_Info.txt", "Gene_info.txt",
                          "Summary_Counts.default_predictions.txt"))
    ]]
  ))
  if (length(missing_columns) > 0L) {
    stop(
      table_name, " is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

mouse_family <- mirna_family[mirna_family[["Species ID"]] == 10090, , drop = FALSE]
mouse_genes <- gene_info[gene_info[["Species ID"]] == 10090, , drop = FALSE]
mouse_predictions <- predictions[
  predictions[["Species ID"]] == 10090,
  ,
  drop = FALSE
]

family_lookup <- mouse_family[, c("miR family", "MiRBase ID"), drop = FALSE]
names(family_lookup) <- c("mirna_family", "mirbase_id")
family_lookup <- unique(family_lookup)

gene_lookup <- mouse_genes[, c("Transcript ID", "Gene ID", "Gene symbol"), drop = FALSE]
names(gene_lookup) <- c("target_transcript_id", "target_entrez_id", "gene_symbol_info")
gene_lookup <- unique(gene_lookup)

standardized <- merge(
  mouse_predictions,
  family_lookup,
  by.x = "miRNA family",
  by.y = "mirna_family",
  all.x = TRUE,
  sort = FALSE
)
standardized <- merge(
  standardized,
  gene_lookup,
  by.x = "Transcript ID",
  by.y = "target_transcript_id",
  all.x = TRUE,
  sort = FALSE
)

target_symbol <- trimws(as.character(standardized[["Gene Symbol"]]))
target_symbol[target_symbol == ""] <- standardized$gene_symbol_info[target_symbol == ""]

output <- data.frame(
  mirna_id = tolower(trimws(as.character(standardized[["Representative miRNA"]]))),
  target_gene_symbol = target_symbol,
  # TargetScan's "Gene ID" column contains Ensembl gene IDs, not Entrez IDs.
  target_gene_id_original = trimws(as.character(standardized$target_entrez_id)),
  target_ensembl_gene_id = sub(
    "\\.[0-9]+$", "", trimws(as.character(standardized$target_entrez_id))
  ),
  target_transcript_id_original = trimws(as.character(standardized[["Transcript ID"]])),
  target_ensembl_transcript_id = sub(
    "\\.[0-9]+$", "", trimws(as.character(standardized[["Transcript ID"]]))
  ),
  evidence_source = "TargetScanMouse",
  evidence_type = "predicted_conserved_target",
  evidence_detail = NA_character_,
  reference_id = NA_character_,
  targetscan_context_score = standardized[["Total context++ score"]],
  targetscan_weighted_context_score = standardized[["Cumulative weighted context++ score"]],
  targetscan_aggregate_pct = standardized[["Aggregate PCT"]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

output <- unique(output)
output <- output[
  nzchar(output$mirna_id) & nzchar(output$target_gene_symbol),
  ,
  drop = FALSE
]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(output, output_file, row.names = FALSE, na = "")

cat("Wrote ", nrow(output), " standardized TargetScanMouse relationships: ",
    output_file, "\n", sep = "")
