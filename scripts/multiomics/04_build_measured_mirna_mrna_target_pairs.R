#!/usr/bin/env Rscript

# Match standardized mouse miRNA-target evidence to measured miRNAs and
# analysis-ready mRNA genes. Matching is case-insensitive for miRNA and gene
# symbols, and uses Entrez IDs or normalized Ensembl gene IDs when available.
#
# Inputs:
#   results/multiomics/miRNA_target_evidence/combined/
#     mirna_mouse_target_evidence.csv
#   results/mirna/expression/rma_normalized_mirna/annotation/
#     mouse_mature_mirna_annotation_mirna.csv
#   results/mrna/analysis_ready/
#     transcript_cluster_annotation_unique_gene_mapped_mrna.csv
#
# Output:
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_target_pairs.csv
#
# Usage:
#   Rscript scripts/multiomics/04_build_measured_mirna_mrna_target_pairs.R

target_file <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "combined",
  "mirna_mouse_target_evidence.csv"
)
mirna_annotation_file <- file.path(
  "results", "mirna", "expression", "rma_normalized_mirna", "annotation",
  "mouse_mature_mirna_annotation_mirna.csv"
)
mrna_annotation_file <- file.path(
  "results", "mrna", "analysis_ready",
  "transcript_cluster_annotation_unique_gene_mapped_mrna.csv"
)
output_dir <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "analysis_ready"
)
output_file <- file.path(output_dir, "mirna_mrna_target_pairs.csv")

required_files <- c(target_file, mirna_annotation_file, mrna_annotation_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required input file(s): ", paste(missing_files, collapse = ", "),
       call. = FALSE)
}

read_table <- function(path) {
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
           na.strings = c("", "NA", "N/A"))
}

clean_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("", "na", "null")] <- NA_character_
  x
}

strip_version <- function(x) {
  x <- trimws(as.character(x))
  sub("\\.[0-9]+$", "", x)
}

targets <- read_table(target_file)
mirna_annotation <- read_table(mirna_annotation_file)
mrna_annotation <- read_table(mrna_annotation_file)

required_target <- c("mirna_id", "target_gene_symbol", "evidence_source")
required_mirna <- c("ProbeSetName", "Transcript ID(Array Design)")
required_mrna <- c("transcript_cluster_id", "ENTREZID", "SYMBOL", "ENSEMBL")
for (required in list(required_target, required_mirna, required_mrna)) {
  table_name <- if (identical(required, required_target)) "target evidence" else
    if (identical(required, required_mirna)) "miRNA annotation" else "mRNA annotation"
  missing <- setdiff(required, names(
    list(targets, mirna_annotation, mrna_annotation)[[
      match(table_name, c("target evidence", "miRNA annotation", "mRNA annotation"))
    ]]
  ))
  if (length(missing) > 0L) {
    stop(table_name, " is missing column(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
}

mirna_map <- data.frame(
  mirna_id = clean_key(mirna_annotation[["Transcript ID(Array Design)"]]),
  mirna_probeset_id = as.character(mirna_annotation$ProbeSetName),
  stringsAsFactors = FALSE
)
mirna_map <- unique(mirna_map[!is.na(mirna_map$mirna_id), , drop = FALSE])

mrna_map <- data.frame(
  transcript_cluster_id = as.character(mrna_annotation$transcript_cluster_id),
  target_entrez_id = clean_key(mrna_annotation$ENTREZID),
  target_gene_symbol = clean_key(mrna_annotation$SYMBOL),
  target_ensembl_gene_id = clean_key(strip_version(mrna_annotation$ENSEMBL)),
  gene_name = as.character(mrna_annotation$GENENAME),
  gene_type = as.character(mrna_annotation$GENETYPE),
  stringsAsFactors = FALSE
)
mrna_map <- unique(mrna_map)

targets$mirna_id <- clean_key(targets$mirna_id)
targets$target_gene_symbol_key <- clean_key(targets$target_gene_symbol)
targets$target_entrez_id_key <- clean_key(targets$target_entrez_id)
if (!"target_ensembl_gene_id" %in% names(targets)) {
  targets$target_ensembl_gene_id <- NA_character_
}
targets$target_ensembl_gene_id_key <- clean_key(strip_version(targets$target_ensembl_gene_id))

matched <- merge(targets, mirna_map, by = "mirna_id", all = FALSE, sort = FALSE)
matched$target_gene_symbol_key <- clean_key(matched$target_gene_symbol)
matched$target_entrez_id_key <- clean_key(matched$target_entrez_id)
matched$target_ensembl_gene_id_key <- clean_key(strip_version(matched$target_ensembl_gene_id))

by_entrez <- !is.na(matched$target_entrez_id_key) &
  matched$target_entrez_id_key %in% mrna_map$target_entrez_id
by_ensembl <- !is.na(matched$target_ensembl_gene_id_key) &
  matched$target_ensembl_gene_id_key %in% mrna_map$target_ensembl_gene_id
by_symbol <- !is.na(matched$target_gene_symbol_key) &
  matched$target_gene_symbol_key %in% mrna_map$target_gene_symbol

# Use one explicit join key per relationship, prioritizing stable IDs and
# falling back to a normalized gene symbol. This avoids requiring every
# database to provide every identifier type.
matched$target_match_type <- ifelse(by_entrez, "Entrez ID",
  ifelse(by_ensembl, "Ensembl gene ID", ifelse(by_symbol, "gene symbol", NA_character_)))
matched$match_key <- ifelse(matched$target_match_type == "Entrez ID",
  paste0("entrez:", matched$target_entrez_id_key),
  ifelse(matched$target_match_type == "Ensembl gene ID",
    paste0("ensembl:", matched$target_ensembl_gene_id_key),
    ifelse(matched$target_match_type == "gene symbol",
      paste0("symbol:", matched$target_gene_symbol_key), NA_character_)))
matched <- matched[!is.na(matched$match_key), , drop = FALSE]

gene_keys <- rbind(
  data.frame(match_key = paste0("entrez:", mrna_map$target_entrez_id), mrna_map,
             stringsAsFactors = FALSE),
  data.frame(match_key = paste0("ensembl:", mrna_map$target_ensembl_gene_id), mrna_map,
             stringsAsFactors = FALSE),
  data.frame(match_key = paste0("symbol:", mrna_map$target_gene_symbol), mrna_map,
             stringsAsFactors = FALSE)
)
gene_keys <- gene_keys[!is.na(gene_keys$match_key) & nzchar(gene_keys$match_key), , drop = FALSE]
gene_keys <- unique(gene_keys[, c("match_key", "transcript_cluster_id", "gene_name", "gene_type"), drop = FALSE])
matched <- merge(matched, gene_keys, by = "match_key", all = FALSE, sort = FALSE)

drop_columns <- grep("_key$|^match_key$", names(matched), value = TRUE)
matched <- matched[, setdiff(names(matched), drop_columns), drop = FALSE]
matched <- unique(matched)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(matched, output_file, row.names = FALSE, na = "")

cat("Measured miRNAs matched: ", length(unique(matched$mirna_id)), "\n", sep = "")
cat("Matched target pairs: ", nrow(matched), "\n", sep = "")
cat("Wrote: ", output_file, "\n", sep = "")
