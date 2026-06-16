#!/usr/bin/env Rscript

# Combine the existing DE, target-database, and paired-association evidence
# into one transparent summary row per selected miRNA--mRNA pair.
#
# This script does not run a new statistical test and does not calculate an
# arbitrary composite score. It preserves the existing statistics and assigns
# descriptive evidence categories using explicit direction and FDR rules.
#
# Inputs:
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_fdr_candidates.csv
#   results/multiomics/paired_association/
#     mirna_mrna_adjusted_association_results.csv
#   results/multiomics/paired_association/covariate_effects/
#     mirna_mrna_association_model_comparison_wide.csv
#
# Outputs:
#   results/multiomics/integrated_pair_evidence/
#     mirna_mrna_integrated_pair_evidence_summary.csv
#     mirna_mrna_directionally_compatible_unconfirmed_pairs.csv
#     mirna_mrna_integrated_pair_evidence_summary.txt
#
# Usage:
#   Rscript scripts/multiomics/10_build_integrated_pair_evidence_summary.R

candidate_file <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "analysis_ready",
  "mirna_mrna_fdr_candidates.csv"
)
primary_association_file <- file.path(
  "results", "multiomics", "paired_association",
  "mirna_mrna_adjusted_association_results.csv"
)
model_comparison_file <- file.path(
  "results", "multiomics", "paired_association", "covariate_effects",
  "mirna_mrna_association_model_comparison_wide.csv"
)
output_dir <- file.path(
  "results", "multiomics", "integrated_pair_evidence"
)
integrated_output_file <- file.path(
  output_dir, "mirna_mrna_integrated_pair_evidence_summary.csv"
)
compatible_output_file <- file.path(
  output_dir,
  "mirna_mrna_directionally_compatible_unconfirmed_pairs.csv"
)
summary_output_file <- file.path(
  output_dir, "mirna_mrna_integrated_pair_evidence_summary.txt"
)

required_inputs <- c(
  candidate_file,
  primary_association_file,
  model_comparison_file
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

read_table <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA", "N/A")
  )
}

candidates <- read_table(candidate_file)
primary_associations <- read_table(primary_association_file)
model_comparisons <- read_table(model_comparison_file)

pair_key <- function(data) {
  paste(
    data$mirna_probeset_id,
    data$transcript_cluster_id,
    data$contrast,
    sep = "||"
  )
}

required_key_columns <- c(
  "mirna_probeset_id", "transcript_cluster_id", "contrast"
)
for (table_name in c(
  "candidates", "primary_associations", "model_comparisons"
)) {
  data <- get(table_name)
  missing_key_columns <- setdiff(required_key_columns, names(data))
  if (length(missing_key_columns) > 0L) {
    stop(
      table_name,
      " lacks required pair-key column(s): ",
      paste(missing_key_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

required_candidate_columns <- c(
  "mirna_id", "mirna_probeset_id", "transcript_cluster_id",
  "mrna_transcript_cluster_id", "gene_name", "gene_type",
  "gene_entrez_id", "gene_symbol", "gene_ensembl_id",
  "mirna_logFC", "mirna_P.Value", "mirna_adj.P.Val",
  "mrna_logFC", "mrna_P.Value", "mrna_adj.P.Val",
  "target_match_type_de", "contrast",
  "mirna_fdr_lt_0_10", "mrna_fdr_lt_0_10",
  "both_fdr_lt_0_10", "inverse_logFC_direction",
  "evidence_record_count", "mirtarbase_record_count",
  "targetscanmouse_record_count", "mirtarbase_supported",
  "targetscanmouse_supported", "both_sources_supported",
  "evidence_sources", "mirtarbase_reference_count",
  "mirtarbase_reference_ids", "mirtarbase_evidence_types",
  "mirtarbase_evidence_details",
  "targetscanmouse_context_score_min",
  "targetscanmouse_context_score_max",
  "targetscanmouse_weighted_context_score_min",
  "targetscanmouse_weighted_context_score_max",
  "targetscanmouse_aggregate_pct_min",
  "targetscanmouse_aggregate_pct_max"
)
missing_candidate_columns <- setdiff(
  required_candidate_columns,
  names(candidates)
)
if (length(missing_candidate_columns) > 0L) {
  stop(
    "Candidate table lacks required column(s): ",
    paste(missing_candidate_columns, collapse = ", "),
    call. = FALSE
  )
}

required_primary_columns <- c(
  "n_paired_animals", "unadjusted_pearson_r",
  "unadjusted_P.Value", "adjusted_beta", "adjusted_SE",
  "adjusted_t", "adjusted_residual_df", "adjusted_partial_r",
  "adjusted_ci_lower", "adjusted_ci_upper",
  "adjusted_P.Value", "adjusted_adj.P.Val",
  "association_direction", "negative_adjusted_association",
  "adjusted_fdr_lt_0_05", "adjusted_fdr_lt_0_10",
  "association_rank"
)
missing_primary_columns <- setdiff(
  required_primary_columns,
  names(primary_associations)
)
if (length(missing_primary_columns) > 0L) {
  stop(
    "Primary association table lacks required column(s): ",
    paste(missing_primary_columns, collapse = ", "),
    call. = FALSE
  )
}

required_comparison_columns <- c(
  "unadjusted_all_43_adj.P.Val",
  "unadjusted_all_43_fdr_lt_0_10",
  "unadjusted_all_43_fdr_lt_0_05",
  "treatment_adjusted_all_43_beta",
  "treatment_adjusted_all_43_partial_r",
  "treatment_adjusted_all_43_P.Value",
  "treatment_adjusted_all_43_adj.P.Val",
  "treatment_adjusted_all_43_association_direction",
  "treatment_adjusted_all_43_fdr_lt_0_10",
  "treatment_adjusted_all_43_fdr_lt_0_05",
  "sex_age_adjusted_all_43_beta",
  "sex_age_adjusted_all_43_partial_r",
  "sex_age_adjusted_all_43_P.Value",
  "sex_age_adjusted_all_43_adj.P.Val",
  "sex_age_adjusted_all_43_association_direction",
  "sex_age_adjusted_all_43_fdr_lt_0_10",
  "sex_age_adjusted_all_43_fdr_lt_0_05",
  "fully_adjusted_all_43_beta",
  "fully_adjusted_all_43_partial_r",
  "fully_adjusted_all_43_P.Value",
  "fully_adjusted_all_43_adj.P.Val",
  "fully_adjusted_all_43_association_direction",
  "fully_adjusted_all_43_fdr_lt_0_10",
  "fully_adjusted_all_43_fdr_lt_0_05",
  "fully_adjusted_mcao_27_n_animals",
  "fully_adjusted_mcao_27_beta",
  "fully_adjusted_mcao_27_partial_r",
  "fully_adjusted_mcao_27_P.Value",
  "fully_adjusted_mcao_27_adj.P.Val",
  "fully_adjusted_mcao_27_association_direction",
  "fully_adjusted_mcao_27_fdr_lt_0_10",
  "fully_adjusted_mcao_27_fdr_lt_0_05",
  "lost_raw_significance_after_treatment_adjustment",
  "lost_raw_significance_after_sex_age_adjustment",
  "lost_raw_significance_after_full_adjustment"
)
missing_comparison_columns <- setdiff(
  required_comparison_columns,
  names(model_comparisons)
)
if (length(missing_comparison_columns) > 0L) {
  stop(
    "Model-comparison table lacks required column(s): ",
    paste(missing_comparison_columns, collapse = ", "),
    call. = FALSE
  )
}

if (
  nrow(candidates) != 56L ||
    nrow(primary_associations) != 56L ||
    nrow(model_comparisons) != 56L
) {
  stop(
    "Expected 56 rows in each input table; found ",
    paste(
      c(
        nrow(candidates),
        nrow(primary_associations),
        nrow(model_comparisons)
      ),
      collapse = ", "
    ),
    ".",
    call. = FALSE
  )
}

candidate_key <- pair_key(candidates)
primary_key <- pair_key(primary_associations)
comparison_key <- pair_key(model_comparisons)
if (
  anyDuplicated(candidate_key) ||
    anyDuplicated(primary_key) ||
    anyDuplicated(comparison_key) ||
    !setequal(candidate_key, primary_key) ||
    !setequal(candidate_key, comparison_key)
) {
  stop(
    "The three inputs do not contain the same 56 unique pair keys.",
    call. = FALSE
  )
}

primary_match <- match(candidate_key, primary_key)
comparison_match <- match(candidate_key, comparison_key)
primary_associations <- primary_associations[
  primary_match,
  ,
  drop = FALSE
]
model_comparisons <- model_comparisons[
  comparison_match,
  ,
  drop = FALSE
]

constant_columns <- c(
  "mirna_id", "mirna_probeset_id", "transcript_cluster_id",
  "gene_name", "gene_entrez_id", "gene_symbol", "contrast",
  "mirna_logFC", "mirna_P.Value", "mirna_adj.P.Val",
  "mrna_logFC", "mrna_P.Value", "mrna_adj.P.Val",
  "inverse_logFC_direction", "mirtarbase_supported",
  "targetscanmouse_supported", "both_sources_supported",
  "evidence_sources"
)
for (column in constant_columns) {
  candidate_value <- as.character(candidates[[column]])
  primary_value <- as.character(primary_associations[[column]])
  comparison_value <- as.character(model_comparisons[[column]])
  candidate_value[is.na(candidate_value)] <- "<NA>"
  primary_value[is.na(primary_value)] <- "<NA>"
  comparison_value[is.na(comparison_value)] <- "<NA>"
  if (
    !identical(candidate_value, primary_value) ||
      !identical(candidate_value, comparison_value)
  ) {
    stop(
      column,
      " does not agree across the three aligned inputs.",
      call. = FALSE
    )
  }
}

if (
  !all(candidates$both_fdr_lt_0_10) ||
    sum(candidates$evidence_record_count) != 59L ||
    !all(primary_associations$n_paired_animals == 43L)
) {
  stop(
    "Candidate-selection or evidence-accounting checks failed.",
    call. = FALSE
  )
}

if (
  !isTRUE(all.equal(
    primary_associations$adjusted_beta,
    model_comparisons$fully_adjusted_all_43_beta,
    tolerance = 1e-12
  )) ||
    !isTRUE(all.equal(
      primary_associations$adjusted_P.Value,
      model_comparisons$fully_adjusted_all_43_P.Value,
      tolerance = 1e-12
    )) ||
    !isTRUE(all.equal(
      primary_associations$adjusted_adj.P.Val,
      model_comparisons$fully_adjusted_all_43_adj.P.Val,
      tolerance = 1e-12
    ))
) {
  stop(
    "Primary and model-comparison fully adjusted results disagree.",
    call. = FALSE
  )
}

target_evidence_category <- ifelse(
  candidates$both_sources_supported,
  "miRTarBase_and_TargetScanMouse",
  ifelse(
    candidates$mirtarbase_supported,
    "miRTarBase_only",
    ifelse(
      candidates$targetscanmouse_supported,
      "TargetScanMouse_only",
      "no_target_source"
    )
  )
)
if (any(target_evidence_category == "no_target_source")) {
  stop(
    "A candidate pair lacks miRTarBase and TargetScanMouse support.",
    call. = FALSE
  )
}

negative_full_association <-
  primary_associations$negative_adjusted_association
negative_mcao_association <-
  model_comparisons$fully_adjusted_mcao_27_association_direction ==
    "negative"
full_association_fdr_supported <-
  primary_associations$adjusted_fdr_lt_0_10
mcao_association_fdr_supported <-
  model_comparisons$fully_adjusted_mcao_27_fdr_lt_0_10

directionally_compatible_full <-
  candidates$inverse_logFC_direction &
    negative_full_association
directionally_compatible_mcao <-
  candidates$inverse_logFC_direction &
    negative_mcao_association
fully_supported_regulatory_pattern <-
  directionally_compatible_full &
    full_association_fdr_supported
directionally_compatible_unconfirmed <-
  directionally_compatible_full &
    !full_association_fdr_supported

integrated_evidence_category <- rep(
  NA_character_,
  nrow(candidates)
)
integrated_evidence_category[
  fully_supported_regulatory_pattern
] <- "inverse_de_negative_association_fdr_lt_0.10"
integrated_evidence_category[
  is.na(integrated_evidence_category) &
    directionally_compatible_unconfirmed
] <- "inverse_de_negative_association_not_fdr_supported"
integrated_evidence_category[
  is.na(integrated_evidence_category) &
    candidates$inverse_logFC_direction
] <- "inverse_de_without_negative_adjusted_association"
integrated_evidence_category[
  is.na(integrated_evidence_category) &
    full_association_fdr_supported &
    negative_full_association
] <- "negative_association_fdr_lt_0.10_without_inverse_de"
integrated_evidence_category[
  is.na(integrated_evidence_category) &
    full_association_fdr_supported
] <- "positive_association_fdr_lt_0.10_without_inverse_de"
integrated_evidence_category[
  is.na(integrated_evidence_category)
] <- "same_de_direction_without_adjusted_association_fdr_support"

plain_english_evidence_summary <- ifelse(
  integrated_evidence_category ==
    "inverse_de_negative_association_fdr_lt_0.10",
  paste(
    "Inverse DE directions and a negative adjusted association",
    "that passes pair-level FDR < 0.10."
  ),
  ifelse(
    integrated_evidence_category ==
      "inverse_de_negative_association_not_fdr_supported",
    paste(
      "Inverse DE directions and a negative adjusted association,",
      "but the association does not pass pair-level FDR < 0.10."
    ),
    ifelse(
      integrated_evidence_category ==
        "inverse_de_without_negative_adjusted_association",
      paste(
        "Inverse DE directions, but the fully adjusted association",
        "is not negative."
      ),
      ifelse(
        integrated_evidence_category ==
          "negative_association_fdr_lt_0.10_without_inverse_de",
        paste(
          "Negative adjusted association passes pair-level FDR < 0.10,",
          "but the DE directions are not inverse."
        ),
        ifelse(
          integrated_evidence_category ==
            "positive_association_fdr_lt_0.10_without_inverse_de",
          paste(
            "Positive adjusted association passes pair-level FDR < 0.10,",
            "and the DE directions are not inverse."
          ),
          paste(
            "The DE directions are the same and the adjusted association",
            "does not pass pair-level FDR < 0.10."
          )
        )
      )
    )
  )
)

category_columns <- data.frame(
  integrated_evidence_category = integrated_evidence_category,
  plain_english_evidence_summary = plain_english_evidence_summary,
  target_evidence_category = target_evidence_category,
  directionally_compatible_full = directionally_compatible_full,
  directionally_compatible_mcao = directionally_compatible_mcao,
  directionally_compatible_unconfirmed =
    directionally_compatible_unconfirmed,
  fully_supported_regulatory_pattern =
    fully_supported_regulatory_pattern,
  full_association_fdr_supported =
    full_association_fdr_supported,
  mcao_association_fdr_supported =
    mcao_association_fdr_supported,
  stringsAsFactors = FALSE
)

primary_output_columns <- c(
  "n_paired_animals", "unadjusted_pearson_r",
  "unadjusted_P.Value", "adjusted_beta", "adjusted_SE",
  "adjusted_t", "adjusted_residual_df", "adjusted_partial_r",
  "adjusted_ci_lower", "adjusted_ci_upper",
  "adjusted_P.Value", "adjusted_adj.P.Val",
  "association_direction", "negative_adjusted_association",
  "adjusted_fdr_lt_0_05", "adjusted_fdr_lt_0_10",
  "association_rank"
)
comparison_output_columns <- required_comparison_columns

integrated_results <- cbind(
  category_columns,
  candidates,
  primary_associations[
    ,
    primary_output_columns,
    drop = FALSE
  ],
  model_comparisons[
    ,
    comparison_output_columns,
    drop = FALSE
  ]
)

category_order <- match(
  integrated_results$integrated_evidence_category,
  c(
    "inverse_de_negative_association_fdr_lt_0.10",
    "inverse_de_negative_association_not_fdr_supported",
    "inverse_de_without_negative_adjusted_association",
    "negative_association_fdr_lt_0.10_without_inverse_de",
    "positive_association_fdr_lt_0.10_without_inverse_de",
    "same_de_direction_without_adjusted_association_fdr_support"
  )
)
output_order <- order(
  category_order,
  integrated_results$adjusted_P.Value,
  integrated_results$mirna_id,
  integrated_results$gene_symbol,
  na.last = TRUE
)
integrated_results <- integrated_results[
  output_order,
  ,
  drop = FALSE
]
row.names(integrated_results) <- NULL

compatible_results <- integrated_results[
  integrated_results$directionally_compatible_unconfirmed,
  ,
  drop = FALSE
]
compatible_results <- compatible_results[
  order(
    compatible_results$adjusted_P.Value,
    compatible_results$mirna_id,
    compatible_results$gene_symbol
  ),
  ,
  drop = FALSE
]
row.names(compatible_results) <- NULL

expected_category_counts <- c(
  inverse_de_negative_association_not_fdr_supported = 10L,
  inverse_de_without_negative_adjusted_association = 3L,
  same_de_direction_without_adjusted_association_fdr_support = 43L
)
observed_category_counts <- table(
  integrated_results$integrated_evidence_category
)
category_counts_match <- all(vapply(
  names(expected_category_counts),
  function(category_name) {
    observed_count <- observed_category_counts[category_name]
    length(observed_count) == 1L &&
      !is.na(observed_count) &&
      unname(observed_count) ==
        unname(expected_category_counts[category_name])
  },
  logical(1)
))
if (
  nrow(integrated_results) != 56L ||
    nrow(compatible_results) != 10L ||
    sum(integrated_results$inverse_logFC_direction) != 13L ||
    any(integrated_results$fully_supported_regulatory_pattern) ||
    any(integrated_results$full_association_fdr_supported) ||
    any(integrated_results$mcao_association_fdr_supported) ||
    !category_counts_match ||
    !all(compatible_results$inverse_logFC_direction) ||
    !all(compatible_results$negative_adjusted_association) ||
    any(compatible_results$adjusted_fdr_lt_0_10)
) {
  stop(
    "Integrated evidence categories failed expected validation checks.",
    call. = FALSE
  )
}

source_counts <- table(integrated_results$target_evidence_category)
compatible_source_counts <- table(
  compatible_results$target_evidence_category
)
category_counts <- table(
  integrated_results$integrated_evidence_category
)

summary_lines <- c(
  "Integrated miRNA--mRNA pair evidence summary",
  "",
  "This script organizes existing evidence and runs no new statistical test.",
  "",
  paste0("Candidate input: ", candidate_file),
  paste0("Primary association input: ", primary_association_file),
  paste0("Model-comparison input: ", model_comparison_file),
  paste0("Selected measured pairs: ", nrow(integrated_results)),
  paste0(
    "Unique miRNAs: ",
    length(unique(integrated_results$mirna_id))
  ),
  paste0(
    "Unique genes: ",
    length(unique(integrated_results$gene_entrez_id))
  ),
  paste0(
    "Inverse-DE pairs: ",
    sum(integrated_results$inverse_logFC_direction)
  ),
  paste0(
    "Inverse-DE pairs with a negative fully adjusted slope: ",
    sum(integrated_results$directionally_compatible_full)
  ),
  paste0(
    "Fully adjusted associations with FDR < 0.10: ",
    sum(integrated_results$full_association_fdr_supported)
  ),
  paste0(
    "MCAO-focused adjusted associations with FDR < 0.10: ",
    sum(integrated_results$mcao_association_fdr_supported)
  ),
  paste0(
    "Directionally compatible but association-unconfirmed pairs: ",
    nrow(compatible_results)
  ),
  "",
  "Target-evidence categories:"
)
for (source_name in names(source_counts)) {
  summary_lines <- c(
    summary_lines,
    paste0("- ", source_name, ": ", source_counts[[source_name]])
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "Integrated evidence categories:"
)
for (category_name in names(category_counts)) {
  summary_lines <- c(
    summary_lines,
    paste0("- ", category_name, ": ", category_counts[[category_name]])
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "Target evidence among the 10 directionally compatible pairs:"
)
for (source_name in names(compatible_source_counts)) {
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ", source_name, ": ",
      compatible_source_counts[[source_name]]
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  paste0("Wrote: ", integrated_output_file),
  paste0("Wrote: ", compatible_output_file),
  "",
  paste(
    "Interpretation boundary: directionally compatible pairs are",
    "exploratory and association-unconfirmed. These categories do not",
    "establish direct miRNA regulation, causality, or mediation."
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  integrated_results,
  integrated_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  compatible_results,
  compatible_output_file,
  row.names = FALSE,
  na = ""
)
writeLines(summary_lines, summary_output_file)

cat("Selected measured pairs: ", nrow(integrated_results), "\n", sep = "")
cat(
  "Inverse-DE pairs: ",
  sum(integrated_results$inverse_logFC_direction),
  "\n",
  sep = ""
)
cat(
  "Directionally compatible but association-unconfirmed pairs: ",
  nrow(compatible_results),
  "\n",
  sep = ""
)
cat(
  "Fully adjusted associations with FDR < 0.10: ",
  sum(integrated_results$full_association_fdr_supported),
  "\n",
  sep = ""
)
cat("Wrote: ", integrated_output_file, "\n", sep = "")
cat("Wrote: ", compatible_output_file, "\n", sep = "")
cat("Wrote: ", summary_output_file, "\n", sep = "")
