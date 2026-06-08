#!/usr/bin/env Rscript

# Collapse evidence-level miRNA--mRNA rows to one row per exact measured
# miRNA probeset--mRNA transcript cluster--contrast combination.
#
# The script keeps differential-expression statistics at their measured
# feature level and summarizes miRTarBase and TargetScanMouse evidence in
# separate columns. It does not rank pairs, run association tests, perform
# pathway analysis, or make biological interpretations.
#
# Input:
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_de_pair_evidence_summary.csv
#
# Output:
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_unique_pair_summary.csv
#     mirna_mrna_fdr_candidates.csv
#
# Usage:
#   Rscript scripts/multiomics/07_collapse_target_evidence_to_unique_pairs.R

input_file <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "analysis_ready",
  "mirna_mrna_de_pair_evidence_summary.csv"
)
output_dir <- dirname(input_file)
output_file <- file.path(
  output_dir, "mirna_mrna_unique_pair_summary.csv"
)
candidate_output_file <- file.path(
  output_dir, "mirna_mrna_fdr_candidates.csv"
)

if (!file.exists(input_file)) {
  stop("Missing evidence-level pair table: ", input_file, call. = FALSE)
}

read_table <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "N/A")
  )
}

clean_key <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("", "na", "null")] <- NA_character_
  x
}

is_present <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

assert_constant_by_group <- function(data, group_id, columns) {
  for (column in columns) {
    value <- as.character(data[[column]])
    value[is.na(value)] <- "<NA>"
    distinct_record <- !duplicated(paste(group_id, value, sep = "\r"))
    distinct_count <- tabulate(
      group_id[distinct_record],
      nbins = max(group_id)
    )
    varying_groups <- sum(distinct_count > 1L)
    if (varying_groups > 0L) {
      stop(
        column, " varies within ", varying_groups,
        " exact measured pair--contrast group(s).",
        call. = FALSE
      )
    }
  }
}

collapse_group_values <- function(values, group_id, n_groups,
                                  separator = " || ") {
  values <- trimws(as.character(values))
  keep <- is_present(values)
  output <- rep(NA_character_, n_groups)
  if (!any(keep)) {
    return(output)
  }

  kept_groups <- group_id[keep]
  kept_values <- values[keep]
  distinct_record <- !duplicated(
    paste(kept_groups, kept_values, sep = "\r")
  )
  kept_groups <- kept_groups[distinct_record]
  kept_values <- kept_values[distinct_record]

  distinct_count <- tabulate(kept_groups, nbins = n_groups)
  one_value_groups <- which(distinct_count == 1L)
  first_position <- match(one_value_groups, kept_groups)
  output[one_value_groups] <- kept_values[first_position]

  multiple_value_groups <- which(distinct_count > 1L)
  if (length(multiple_value_groups) > 0L) {
    use_multiple <- kept_groups %in% multiple_value_groups
    grouped_values <- split(
      kept_values[use_multiple],
      kept_groups[use_multiple],
      drop = TRUE
    )
    collapsed <- vapply(
      grouped_values,
      function(x) paste(sort(unique(x)), collapse = separator),
      character(1)
    )
    output[as.integer(names(collapsed))] <- collapsed
  }

  output
}

count_unique_group_values <- function(values, group_id, n_groups) {
  values <- trimws(as.character(values))
  keep <- is_present(values)
  if (!any(keep)) {
    return(integer(n_groups))
  }
  distinct_record <- !duplicated(
    paste(group_id[keep], values[keep], sep = "\r")
  )
  tabulate(
    group_id[keep][distinct_record],
    nbins = n_groups
  )
}

group_numeric_extreme <- function(values, group_id, n_groups,
                                  method = c("min", "max")) {
  method <- match.arg(method)
  values <- suppressWarnings(as.numeric(values))
  keep <- is.finite(values)
  output <- rep(NA_real_, n_groups)
  if (!any(keep)) {
    return(output)
  }

  kept_groups <- group_id[keep]
  kept_values <- values[keep]
  group_count <- tabulate(kept_groups, nbins = n_groups)
  first_record <- !duplicated(kept_groups)
  output[kept_groups[first_record]] <- kept_values[first_record]

  repeated_groups <- which(group_count > 1L)
  if (length(repeated_groups) > 0L) {
    use_repeated <- kept_groups %in% repeated_groups
    grouped_values <- split(
      kept_values[use_repeated],
      kept_groups[use_repeated],
      drop = TRUE
    )
    extreme_function <- if (method == "min") min else max
    extremes <- vapply(grouped_values, extreme_function, numeric(1))
    output[as.integer(names(extremes))] <- extremes
  }

  output
}

pairs <- read_table(input_file)

key_columns <- c(
  "mirna_probeset_id", "transcript_cluster_id", "contrast"
)
constant_columns <- c(
  "mirna_id", "mirna_probeset_id",
  "transcript_cluster_id", "mrna_transcript_cluster_id",
  "gene_name", "gene_type",
  "mrna_ENTREZID", "mrna_SYMBOL", "mrna_ENSEMBL",
  "mirna_logFC", "mirna_P.Value", "mirna_adj.P.Val",
  "mrna_logFC", "mrna_P.Value", "mrna_adj.P.Val",
  "target_match_type_de", "contrast",
  "mirna_fdr_lt_0_10", "mrna_fdr_lt_0_10",
  "both_fdr_lt_0_10", "inverse_logFC_direction",
  "candidate_category"
)
evidence_columns <- c(
  "evidence_source", "evidence_type", "evidence_detail", "reference_id",
  "target_gene_symbol", "target_entrez_id",
  "target_gene_id_original", "target_ensembl_gene_id",
  "target_transcript_id_original", "target_ensembl_transcript_id",
  "targetscan_context_score", "targetscan_weighted_context_score",
  "targetscan_aggregate_pct", "target_match_type"
)
required_columns <- unique(c(
  key_columns, constant_columns, evidence_columns
))
missing_columns <- setdiff(required_columns, names(pairs))
if (length(missing_columns) > 0L) {
  stop(
    "Input is missing column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(pairs) == 0L) {
  stop("Input evidence-level pair table is empty.", call. = FALSE)
}

supported_sources <- c("miRTarBase", "TargetScanMouse")
unknown_sources <- setdiff(unique(pairs$evidence_source), supported_sources)
if (length(unknown_sources) > 0L) {
  stop(
    "Input contains unsupported evidence_source values: ",
    paste(unknown_sources, collapse = ", "),
    call. = FALSE
  )
}

exact_cluster_match <- clean_key(pairs$transcript_cluster_id) ==
  clean_key(pairs$mrna_transcript_cluster_id)
if (any(is.na(exact_cluster_match)) || !all(exact_cluster_match)) {
  stop(
    "Input contains non-identical target-pair and DE transcript-cluster IDs.",
    call. = FALSE
  )
}

pair_key <- do.call(
  paste,
  c(lapply(pairs[key_columns], clean_key), sep = "||")
)
if (any(is.na(pair_key) | !nzchar(pair_key))) {
  stop("Cannot create a complete measured pair--contrast key.", call. = FALSE)
}

unique_keys <- unique(pair_key)
group_id <- match(pair_key, unique_keys)
n_groups <- length(unique_keys)
assert_constant_by_group(pairs, group_id, constant_columns)

first_record <- !duplicated(group_id)
unique_pairs <- pairs[first_record, constant_columns, drop = FALSE]
row.names(unique_pairs) <- NULL
names(unique_pairs)[names(unique_pairs) == "mrna_ENTREZID"] <- "gene_entrez_id"
names(unique_pairs)[names(unique_pairs) == "mrna_SYMBOL"] <- "gene_symbol"
names(unique_pairs)[names(unique_pairs) == "mrna_ENSEMBL"] <- "gene_ensembl_id"

mirna_rows <- pairs$evidence_source == "miRTarBase"
targetscan_rows <- pairs$evidence_source == "TargetScanMouse"
mirna_group_id <- group_id[mirna_rows]
targetscan_group_id <- group_id[targetscan_rows]

unique_pairs$evidence_record_count <- tabulate(group_id, nbins = n_groups)
unique_pairs$mirtarbase_record_count <- tabulate(
  mirna_group_id, nbins = n_groups
)
unique_pairs$targetscanmouse_record_count <- tabulate(
  targetscan_group_id, nbins = n_groups
)
unique_pairs$mirtarbase_supported <-
  unique_pairs$mirtarbase_record_count > 0L
unique_pairs$targetscanmouse_supported <-
  unique_pairs$targetscanmouse_record_count > 0L
unique_pairs$both_sources_supported <-
  unique_pairs$mirtarbase_supported &
    unique_pairs$targetscanmouse_supported
unique_pairs$evidence_sources <- collapse_group_values(
  pairs$evidence_source, group_id, n_groups
)

unique_pairs$mirtarbase_reference_count <- count_unique_group_values(
  pairs$reference_id[mirna_rows], mirna_group_id, n_groups
)
unique_pairs$mirtarbase_reference_ids <- collapse_group_values(
  pairs$reference_id[mirna_rows], mirna_group_id, n_groups
)
unique_pairs$mirtarbase_evidence_types <- collapse_group_values(
  pairs$evidence_type[mirna_rows], mirna_group_id, n_groups
)
unique_pairs$mirtarbase_evidence_details <- collapse_group_values(
  pairs$evidence_detail[mirna_rows], mirna_group_id, n_groups
)
unique_pairs$mirtarbase_target_gene_symbols <- collapse_group_values(
  pairs$target_gene_symbol[mirna_rows], mirna_group_id, n_groups
)
unique_pairs$mirtarbase_target_entrez_ids <- collapse_group_values(
  pairs$target_entrez_id[mirna_rows], mirna_group_id, n_groups
)
unique_pairs$mirtarbase_target_match_types <- collapse_group_values(
  pairs$target_match_type[mirna_rows], mirna_group_id, n_groups
)

unique_pairs$targetscanmouse_target_gene_symbols <- collapse_group_values(
  pairs$target_gene_symbol[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_original_gene_ids <- collapse_group_values(
  pairs$target_gene_id_original[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_ensembl_gene_ids <- collapse_group_values(
  pairs$target_ensembl_gene_id[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_original_transcript_ids <- collapse_group_values(
  pairs$target_transcript_id_original[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_ensembl_transcript_ids <- collapse_group_values(
  pairs$target_ensembl_transcript_id[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_transcript_count <- count_unique_group_values(
  pairs$target_ensembl_transcript_id[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_evidence_types <- collapse_group_values(
  pairs$evidence_type[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_target_match_types <- collapse_group_values(
  pairs$target_match_type[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_context_scores <- collapse_group_values(
  pairs$targetscan_context_score[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_context_score_min <- group_numeric_extreme(
  pairs$targetscan_context_score[targetscan_rows],
  targetscan_group_id,
  n_groups,
  method = "min"
)
unique_pairs$targetscanmouse_context_score_max <- group_numeric_extreme(
  pairs$targetscan_context_score[targetscan_rows],
  targetscan_group_id,
  n_groups,
  method = "max"
)
unique_pairs$targetscanmouse_weighted_context_scores <- collapse_group_values(
  pairs$targetscan_weighted_context_score[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_weighted_context_score_min <-
  group_numeric_extreme(
    pairs$targetscan_weighted_context_score[targetscan_rows],
    targetscan_group_id,
    n_groups,
    method = "min"
  )
unique_pairs$targetscanmouse_weighted_context_score_max <-
  group_numeric_extreme(
    pairs$targetscan_weighted_context_score[targetscan_rows],
    targetscan_group_id,
    n_groups,
    method = "max"
  )
unique_pairs$targetscanmouse_aggregate_pcts <- collapse_group_values(
  pairs$targetscan_aggregate_pct[targetscan_rows],
  targetscan_group_id,
  n_groups
)
unique_pairs$targetscanmouse_aggregate_pct_min <- group_numeric_extreme(
  pairs$targetscan_aggregate_pct[targetscan_rows],
  targetscan_group_id,
  n_groups,
  method = "min"
)
unique_pairs$targetscanmouse_aggregate_pct_max <- group_numeric_extreme(
  pairs$targetscan_aggregate_pct[targetscan_rows],
  targetscan_group_id,
  n_groups,
  method = "max"
)

if (anyDuplicated(unique_pairs[key_columns])) {
  stop("Output still contains duplicate measured pair--contrast keys.",
       call. = FALSE)
}
if (nrow(unique_pairs) != n_groups) {
  stop("Output row count does not equal the unique input key count.",
       call. = FALSE)
}
if (!all(
  unique_pairs$evidence_record_count ==
    unique_pairs$mirtarbase_record_count +
      unique_pairs$targetscanmouse_record_count
)) {
  stop("Source-specific evidence counts do not sum to total evidence counts.",
       call. = FALSE)
}
if (sum(unique_pairs$mirtarbase_record_count) != sum(mirna_rows)) {
  stop("miRTarBase evidence-row counts were not preserved.", call. = FALSE)
}
if (sum(unique_pairs$targetscanmouse_record_count) != sum(targetscan_rows)) {
  stop("TargetScanMouse evidence-row counts were not preserved.",
       call. = FALSE)
}

contrast_order <- match(unique_pairs$contrast, c(
  "MCAO1hr_vs_Sham",
  "MCAO3hr_vs_Sham",
  "MCAO3hr_vs_MCAO1hr"
))
output_order <- order(
  contrast_order,
  unique_pairs$mirna_id,
  unique_pairs$gene_symbol,
  unique_pairs$mirna_probeset_id,
  unique_pairs$transcript_cluster_id,
  na.last = TRUE
)
unique_pairs <- unique_pairs[output_order, , drop = FALSE]
row.names(unique_pairs) <- NULL

fdr_candidates <- unique_pairs[
  unique_pairs$mirna_fdr_lt_0_10 &
    unique_pairs$mrna_fdr_lt_0_10,
  ,
  drop = FALSE
]
row.names(fdr_candidates) <- NULL

if (anyDuplicated(fdr_candidates[key_columns])) {
  stop("FDR candidate output contains duplicate measured pair--contrast keys.",
       call. = FALSE)
}
if (nrow(fdr_candidates) > 0L &&
    (!all(fdr_candidates$both_fdr_lt_0_10) ||
     !all(fdr_candidates$candidate_category ==
          "both_layers_fdr_lt_0.10"))) {
  stop("FDR candidate output contains a row that fails the selection rule.",
       call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(unique_pairs, output_file, row.names = FALSE, na = "")
write.csv(fdr_candidates, candidate_output_file, row.names = FALSE, na = "")

cat("Input evidence rows: ", nrow(pairs), "\n", sep = "")
cat("Unique measured pair--contrast rows: ", nrow(unique_pairs), "\n", sep = "")
cat(
  "Collapsed evidence rows: ",
  nrow(pairs) - nrow(unique_pairs),
  "\n",
  sep = ""
)
cat("Wrote: ", output_file, "\n", sep = "")
cat(
  "FDR candidates with both layers below 0.10: ",
  nrow(fdr_candidates),
  "\n",
  sep = ""
)
cat("Wrote: ", candidate_output_file, "\n", sep = "")
cat("\nEvidence-source support by unique pair--contrast:\n")
print(table(
  miRTarBase = unique_pairs$mirtarbase_supported,
  TargetScanMouse = unique_pairs$targetscanmouse_supported
))
cat("\nCandidate categories by contrast:\n")
print(table(unique_pairs$contrast, unique_pairs$candidate_category))
