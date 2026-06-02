#!/usr/bin/env Rscript

# Join target-supported miRNA--mRNA relationships to the three treatment
# differential-expression contrasts. The output is long format: a measured
# miRNA--transcript-cluster pair can appear once per evidence record in each
# contrast. Database evidence and one-to-many target mappings are preserved.
#
# Inputs:
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_target_pairs.csv
#   results/mirna/differential_expression/rma_normalized_mirna/
#     treatment_sex_age_limma/dabg20_*_mirna.csv
#   results/mrna/differential_expression/treatment_sex_age_limma/
#     transcript_cluster_*_mrna.csv
#
# Output:
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_de_supported_target_pairs.csv
#
# Usage:
#   Rscript scripts/multiomics/05_join_target_pairs_to_de_results.R

target_file <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "analysis_ready",
  "mirna_mrna_target_pairs.csv"
)
mirna_dir <- file.path(
  "results", "mirna", "differential_expression", "rma_normalized_mirna",
  "treatment_sex_age_limma"
)
mrna_dir <- file.path(
  "results", "mrna", "differential_expression", "treatment_sex_age_limma"
)
output_dir <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "analysis_ready"
)
output_file <- file.path(
  output_dir, "mirna_mrna_de_supported_target_pairs.csv"
)

contrasts <- c(
  "MCAO1hr_vs_Sham",
  "MCAO3hr_vs_Sham",
  "MCAO3hr_vs_MCAO1hr"
)

mirna_files <- setNames(
  file.path(mirna_dir, paste0("dabg20_", contrasts, "_mirna.csv")),
  contrasts
)
mrna_files <- setNames(
  file.path(mrna_dir, paste0("transcript_cluster_", contrasts, "_mrna.csv")),
  contrasts
)

required_files <- c(target_file, mirna_files, mrna_files)
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

target_pairs <- read_table(target_file)
required_target <- c(
  "mirna_id", "mirna_probeset_id", "target_gene_symbol",
  "target_entrez_id", "target_ensembl_gene_id", "evidence_source",
  "target_match_type", "transcript_cluster_id"
)
missing_target <- setdiff(required_target, names(target_pairs))
if (length(missing_target) > 0L) {
  stop("Target-pair table is missing column(s): ",
       paste(missing_target, collapse = ", "), call. = FALSE)
}

target_pairs$mirna_id_key <- clean_key(target_pairs$mirna_id)
target_pairs$mirna_probeset_id_key <- clean_key(target_pairs$mirna_probeset_id)
target_pairs$target_entrez_id_key <- clean_key(target_pairs$target_entrez_id)
target_pairs$target_ensembl_gene_id_key <- clean_key(
  strip_version(target_pairs$target_ensembl_gene_id)
)
target_pairs$target_gene_symbol_key <- clean_key(target_pairs$target_gene_symbol)
target_pairs$transcript_cluster_id_key <- clean_key(
  target_pairs$transcript_cluster_id
)

join_one_contrast <- function(contrast) {
  mirna <- read_table(mirna_files[[contrast]])
  mrna <- read_table(mrna_files[[contrast]])

  required_mirna <- c(
    "ProbeSetName", "Transcript ID(Array Design)",
    "logFC", "P.Value", "adj.P.Val"
  )
  required_mrna <- c(
    "transcript_cluster_id", "ENTREZID", "SYMBOL", "ENSEMBL",
    "logFC", "P.Value", "adj.P.Val"
  )
  missing_mirna <- setdiff(required_mirna, names(mirna))
  missing_mrna <- setdiff(required_mrna, names(mrna))
  if (length(missing_mirna) > 0L) {
    stop(contrast, " miRNA DE file is missing column(s): ",
         paste(missing_mirna, collapse = ", "), call. = FALSE)
  }
  if (length(missing_mrna) > 0L) {
    stop(contrast, " mRNA DE file is missing column(s): ",
         paste(missing_mrna, collapse = ", "), call. = FALSE)
  }

  mirna_de <- data.frame(
    mirna_id_key = clean_key(mirna[["Transcript ID(Array Design)"]]),
    mirna_probeset_id_key = clean_key(mirna$ProbeSetName),
    mirna_logFC = mirna$logFC,
    mirna_P.Value = mirna$P.Value,
    mirna_adj.P.Val = mirna$adj.P.Val,
    stringsAsFactors = FALSE
  )
  mirna_de <- unique(mirna_de)
  mirna_de_key <- paste(
    mirna_de$mirna_id_key,
    mirna_de$mirna_probeset_id_key,
    sep = "||"
  )
  if (anyDuplicated(mirna_de_key)) {
    stop(contrast, " miRNA DE table has duplicate feature identifiers.",
         call. = FALSE)
  }

  # Join the target-pair table back to the exact transcript cluster already
  # selected in script 04. Gene identifiers are retained for validation only;
  # joining by gene can cross-match distinct transcript clusters for genes
  # represented by more than one array feature.
  mrna_de <- data.frame(
    transcript_cluster_id_key = clean_key(mrna$transcript_cluster_id),
    mrna_transcript_cluster_id = as.character(mrna$transcript_cluster_id),
    mrna_ENTREZID = as.character(mrna$ENTREZID),
    mrna_SYMBOL = as.character(mrna$SYMBOL),
    mrna_ENSEMBL = as.character(mrna$ENSEMBL),
    mrna_logFC = mrna$logFC,
    mrna_P.Value = mrna$P.Value,
    mrna_adj.P.Val = mrna$adj.P.Val,
    stringsAsFactors = FALSE
  )
  if (any(is.na(mrna_de$transcript_cluster_id_key))) {
    stop(contrast, " mRNA DE table has missing transcript-cluster IDs.",
         call. = FALSE)
  }
  if (anyDuplicated(mrna_de$transcript_cluster_id_key)) {
    stop(contrast, " mRNA DE table has duplicate transcript-cluster IDs.",
         call. = FALSE)
  }
  mrna_de$target_match_type_de <- "transcript cluster ID"

  pairs <- merge(
    target_pairs,
    mirna_de,
    by = c("mirna_id_key", "mirna_probeset_id_key"),
    all = FALSE,
    sort = FALSE
  )

  missing_clusters <- setdiff(
    unique(pairs$transcript_cluster_id_key),
    mrna_de$transcript_cluster_id_key
  )
  if (length(missing_clusters) > 0L) {
    stop(
      contrast, " target pairs contain transcript clusters absent from ",
      "the mRNA DE table: ", paste(head(missing_clusters, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  pairs <- merge(
    pairs,
    mrna_de,
    by = "transcript_cluster_id_key",
    all = FALSE,
    sort = FALSE,
    suffixes = c("", "_mRNA")
  )

  exact_cluster_match <- clean_key(pairs$transcript_cluster_id) ==
    clean_key(pairs$mrna_transcript_cluster_id)
  if (any(is.na(exact_cluster_match)) || !all(exact_cluster_match)) {
    stop(contrast, " produced a non-identical transcript-cluster join.",
         call. = FALSE)
  }

  # Validate the gene identifier that originally connected the external target
  # record to this measured transcript cluster in script 04.
  match_types <- c("Entrez ID", "Ensembl gene ID", "gene symbol")
  unknown_match_type <- !pairs$target_match_type %in% match_types
  if (any(unknown_match_type)) {
    stop(
      contrast, " target pairs contain unsupported target_match_type values: ",
      paste(unique(pairs$target_match_type[unknown_match_type]), collapse = ", "),
      call. = FALSE
    )
  }

  gene_identifier_match <- rep(FALSE, nrow(pairs))
  by_entrez <- pairs$target_match_type == "Entrez ID"
  by_ensembl <- pairs$target_match_type == "Ensembl gene ID"
  by_symbol <- pairs$target_match_type == "gene symbol"
  gene_identifier_match[by_entrez] <-
    pairs$target_entrez_id_key[by_entrez] ==
      clean_key(pairs$mrna_ENTREZID[by_entrez])
  gene_identifier_match[by_ensembl] <-
    pairs$target_ensembl_gene_id_key[by_ensembl] ==
      clean_key(strip_version(pairs$mrna_ENSEMBL[by_ensembl]))
  gene_identifier_match[by_symbol] <-
    pairs$target_gene_symbol_key[by_symbol] ==
      clean_key(pairs$mrna_SYMBOL[by_symbol])
  if (any(is.na(gene_identifier_match)) || !all(gene_identifier_match)) {
    stop(
      contrast, " failed gene-identifier validation after the exact ",
      "transcript-cluster join.", call. = FALSE
    )
  }

  pairs$contrast <- contrast
  pairs
}

matched_by_contrast <- lapply(contrasts, join_one_contrast)
matched <- do.call(rbind, matched_by_contrast)

drop_columns <- grep("_key$", names(matched), value = TRUE)
matched <- matched[, setdiff(names(matched), drop_columns), drop = FALSE]
matched <- unique(matched)

if (!all(clean_key(matched$transcript_cluster_id) ==
         clean_key(matched$mrna_transcript_cluster_id))) {
  stop("Final output contains non-identical transcript-cluster IDs.",
       call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(matched, output_file, row.names = FALSE, na = "")

cat("Wrote ", nrow(matched), " DE-supported pair rows across ",
    length(contrasts), " contrasts: ", output_file, "\n", sep = "")
print(table(matched$contrast))
cat("Validated exact transcript-cluster joins: ", nrow(matched), "\n", sep = "")
