#!/usr/bin/env Rscript

# Build a validated biological annotation table for the Clariom S Mouse HT
# transcript clusters in the normalized mRNA expression matrix.
#
# The output uses long format: a transcript cluster with multiple gene mappings
# receives one row per mapping. No mapping is silently reduced to the first
# result. Transcript clusters without a gene assignment are retained with an
# explicit mapping status.
#
# Prerequisite:
#   Rscript scripts/mrna/01_array_qc_normalize_mrna.R
#
# Required packages:
#   AnnotationDbi
#   clariomsmousehttranscriptcluster.db
#
# Outputs:
#   results/mrna/annotation/
#     transcript_cluster_gene_annotation_mrna_150002.csv
#     unmatched_transcript_clusters_mrna_150002.csv
#     transcript_cluster_annotation_summary_mrna_150002.txt
#
# Run from the repository root:
#   Rscript scripts/mrna/02_annotate_transcript_clusters_mrna.R

expression_file <- file.path(
  "results", "mrna", "normalized",
  "rma_expression_mrna_150002.csv"
)
output_dir <- file.path("results", "mrna", "annotation")
annotation_output_file <- file.path(
  output_dir,
  "transcript_cluster_gene_annotation_mrna_150002.csv"
)
unmatched_output_file <- file.path(
  output_dir,
  "unmatched_transcript_clusters_mrna_150002.csv"
)
summary_output_file <- file.path(
  output_dir,
  "transcript_cluster_annotation_summary_mrna_150002.txt"
)

required_packages <- c(
  "AnnotationDbi",
  "clariomsmousehttranscriptcluster.db"
)
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}
if (!file.exists(expression_file)) {
  stop(
    "Normalized mRNA expression matrix not found: ",
    expression_file,
    call. = FALSE
  )
}

expression <- read.csv(
  expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!"transcript_cluster_id" %in% names(expression)) {
  stop(
    "Expression matrix lacks transcript_cluster_id.",
    call. = FALSE
  )
}

transcript_cluster_ids <- as.character(expression$transcript_cluster_id)
if (anyNA(transcript_cluster_ids) || any(!nzchar(transcript_cluster_ids))) {
  stop("Expression matrix contains missing transcript-cluster IDs.", call. = FALSE)
}
if (anyDuplicated(transcript_cluster_ids)) {
  stop("Expression matrix contains duplicated transcript-cluster IDs.", call. = FALSE)
}

annotation_db <-
  clariomsmousehttranscriptcluster.db::clariomsmousehttranscriptcluster.db
database_probe_ids <- AnnotationDbi::keys(
  annotation_db,
  keytype = "PROBEID"
)
matched_probe_ids <- intersect(
  transcript_cluster_ids,
  database_probe_ids
)
unmatched_probe_ids <- setdiff(
  transcript_cluster_ids,
  database_probe_ids
)

gene_annotation <- AnnotationDbi::select(
  annotation_db,
  keys = matched_probe_ids,
  keytype = "PROBEID",
  columns = c(
    "ENTREZID",
    "SYMBOL",
    "GENENAME",
    "GENETYPE",
    "ENSEMBL"
  )
)
gene_annotation <- unique(as.data.frame(
  gene_annotation,
  stringsAsFactors = FALSE,
  check.names = FALSE
))
names(gene_annotation)[names(gene_annotation) == "PROBEID"] <-
  "transcript_cluster_id"

gene_annotation$mapping_status <- ifelse(
  is.na(gene_annotation$ENTREZID) &
    is.na(gene_annotation$SYMBOL),
  "database_without_gene",
  "mapped_to_gene"
)

unmatched_annotation <- data.frame(
  transcript_cluster_id = unmatched_probe_ids,
  ENTREZID = NA_character_,
  SYMBOL = NA_character_,
  GENENAME = NA_character_,
  GENETYPE = NA_character_,
  ENSEMBL = NA_character_,
  mapping_status = "not_in_annotation_database",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

annotation <- rbind(gene_annotation, unmatched_annotation)
annotation$expression_row <- match(
  annotation$transcript_cluster_id,
  transcript_cluster_ids
)

mapping_key <- ifelse(
  !is.na(annotation$ENTREZID) & nzchar(annotation$ENTREZID),
  paste0("ENTREZ:", annotation$ENTREZID),
  ifelse(
    !is.na(annotation$SYMBOL) & nzchar(annotation$SYMBOL),
    paste0("SYMBOL:", annotation$SYMBOL),
    NA_character_
  )
)
mapping_count_by_probe <- tapply(
  mapping_key,
  annotation$transcript_cluster_id,
  function(values) length(unique(values[!is.na(values)]))
)
annotation$gene_mapping_count <- as.integer(
  mapping_count_by_probe[annotation$transcript_cluster_id]
)
annotation$is_one_to_many <- annotation$gene_mapping_count > 1L

annotation <- annotation[
  order(
    annotation$expression_row,
    annotation$ENTREZID,
    annotation$SYMBOL,
    na.last = TRUE
  ),
  c(
    "expression_row",
    "transcript_cluster_id",
    "mapping_status",
    "gene_mapping_count",
    "is_one_to_many",
    "ENTREZID",
    "SYMBOL",
    "GENENAME",
    "GENETYPE",
    "ENSEMBL"
  ),
  drop = FALSE
]
rownames(annotation) <- NULL

represented_ids <- unique(annotation$transcript_cluster_id)
if (!setequal(represented_ids, transcript_cluster_ids)) {
  stop("Not every expression transcript cluster is represented.", call. = FALSE)
}
if (anyNA(annotation$expression_row)) {
  stop("Annotation contains an ID absent from the expression matrix.", call. = FALSE)
}
if (anyDuplicated(annotation)) {
  stop("Annotation table contains exact duplicate rows.", call. = FALSE)
}

unmatched_table <- data.frame(
  expression_row = match(unmatched_probe_ids, transcript_cluster_ids),
  transcript_cluster_id = unmatched_probe_ids,
  is_affx_control = grepl(
    "^AFFX", unmatched_probe_ids, ignore.case = TRUE
  ),
  stringsAsFactors = FALSE
)
unmatched_table <- unmatched_table[
  order(unmatched_table$expression_row),
  ,
  drop = FALSE
]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  annotation,
  annotation_output_file,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
write.csv(
  unmatched_table,
  unmatched_output_file,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

status_counts <- table(annotation$mapping_status)
one_to_many_ids <- unique(
  annotation$transcript_cluster_id[annotation$is_one_to_many]
)
writeLines(
  c(
    "Dataset: Clariom S Mouse HT mRNA 150002",
    paste("Expression input:", expression_file),
    paste(
      "Annotation package: clariomsmousehttranscriptcluster.db",
      as.character(utils::packageVersion(
        "clariomsmousehttranscriptcluster.db"
      ))
    ),
    paste("Expression transcript clusters:", length(transcript_cluster_ids)),
    paste("Matched to database PROBEID:", length(matched_probe_ids)),
    paste("Not in annotation database:", length(unmatched_probe_ids)),
    paste("Output annotation rows:", nrow(annotation)),
    paste("Transcript clusters with one-to-many gene mappings:", length(one_to_many_ids)),
    "",
    "Mapping status rows:",
    paste(
      paste(names(status_counts), as.integer(status_counts), sep = ": "),
      collapse = "\n"
    ),
    "",
    "Validation:",
    "Every expression transcript_cluster_id is represented: TRUE",
    "Duplicate expression transcript_cluster_id values: 0",
    "Exact duplicate annotation rows: 0",
    "One-to-many mappings were preserved in long format."
  ),
  summary_output_file
)

message("Wrote annotation: ", annotation_output_file)
message("Wrote unmatched-ID report: ", unmatched_output_file)
message("Wrote summary: ", summary_output_file)
