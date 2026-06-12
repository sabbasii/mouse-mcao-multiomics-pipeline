#!/usr/bin/env Rscript

# Test adjusted expression associations for the DE-supported miRNA--mRNA
# candidate pairs in the same 43 animals.
#
# For each exact measured miRNA--mRNA pair, the primary model is:
#
#   mRNA expression ~ miRNA expression + treatment + sex + age
#
# The miRNA coefficient tests whether animals with higher miRNA expression
# tend to have higher or lower target-mRNA expression after accounting for
# treatment, sex, and age. Benjamini--Hochberg correction is applied across
# all tested pairs.
#
# Inputs:
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_fdr_candidates.csv
#   results/multiomics/sample_manifest/
#     paired_manifest_mirna_mrna.csv
#   results/mirna/expression/rma_normalized_mirna/annotation/
#     mouse_mature_mirna_expression_mirna.csv
#   results/mrna/analysis_ready/
#     expression_matrix_unique_gene_mapped_mrna.csv
#
# Outputs:
#   results/multiomics/paired_association/
#     mirna_mrna_adjusted_association_results.csv
#     mirna_mrna_adjusted_association_summary.txt
#
# Usage:
#   Rscript scripts/multiomics/08_test_adjusted_paired_associations.R

candidate_file <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "analysis_ready",
  "mirna_mrna_fdr_candidates.csv"
)
manifest_file <- file.path(
  "results", "multiomics", "sample_manifest",
  "paired_manifest_mirna_mrna.csv"
)
mirna_expression_file <- file.path(
  "results", "mirna", "expression", "rma_normalized_mirna",
  "annotation", "mouse_mature_mirna_expression_mirna.csv"
)
mrna_expression_file <- file.path(
  "results", "mrna", "analysis_ready",
  "expression_matrix_unique_gene_mapped_mrna.csv"
)
output_dir <- file.path(
  "results", "multiomics", "paired_association"
)
result_file <- file.path(
  output_dir, "mirna_mrna_adjusted_association_results.csv"
)
summary_file <- file.path(
  output_dir, "mirna_mrna_adjusted_association_summary.txt"
)

required_inputs <- c(
  candidate_file,
  manifest_file,
  mirna_expression_file,
  mrna_expression_file
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
manifest <- read_table(manifest_file)
mirna_expression <- read_table(mirna_expression_file)
mrna_expression <- read_table(mrna_expression_file)

required_candidate_columns <- c(
  "mirna_id", "mirna_probeset_id", "transcript_cluster_id",
  "gene_name", "gene_entrez_id", "gene_symbol", "contrast",
  "mirna_logFC", "mirna_P.Value", "mirna_adj.P.Val",
  "mrna_logFC", "mrna_P.Value", "mrna_adj.P.Val",
  "inverse_logFC_direction", "mirna_fdr_lt_0_10",
  "mrna_fdr_lt_0_10", "both_fdr_lt_0_10",
  "mirtarbase_supported", "targetscanmouse_supported",
  "both_sources_supported", "evidence_sources"
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

required_manifest_columns <- c(
  "animal_id", "treatment", "sex", "age_group",
  "mirna_file_name", "mrna_file_name"
)
missing_manifest_columns <- setdiff(
  required_manifest_columns,
  names(manifest)
)
if (length(missing_manifest_columns) > 0L) {
  stop(
    "Paired manifest lacks required column(s): ",
    paste(missing_manifest_columns, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(candidates) == 0L) {
  stop("Candidate table contains no rows to test.", call. = FALSE)
}
if (nrow(manifest) != 43L) {
  stop(
    "Expected 43 paired animals; found ", nrow(manifest), ".",
    call. = FALSE
  )
}
if (
  anyNA(manifest[required_manifest_columns]) ||
    anyDuplicated(manifest$animal_id) ||
    anyDuplicated(manifest$mirna_file_name) ||
    anyDuplicated(manifest$mrna_file_name)
) {
  stop(
    "Paired manifest contains missing or duplicated required identifiers.",
    call. = FALSE
  )
}

pair_key <- paste(
  candidates$mirna_probeset_id,
  candidates$transcript_cluster_id,
  sep = "||"
)
pair_contrast_key <- paste(pair_key, candidates$contrast, sep = "||")
if (anyDuplicated(pair_contrast_key)) {
  stop(
    "Candidate table contains duplicated measured pair--contrast rows.",
    call. = FALSE
  )
}
if (anyDuplicated(pair_key)) {
  stop(
    paste(
      "The same measured pair occurs in more than one candidate row.",
      "Association tests must be unique at the measured-pair level."
    ),
    call. = FALSE
  )
}
if (
  !all(candidates$mirna_fdr_lt_0_10) ||
    !all(candidates$mrna_fdr_lt_0_10) ||
    !all(candidates$both_fdr_lt_0_10)
) {
  stop(
    "Candidate table contains a row that fails the two-layer FDR rule.",
    call. = FALSE
  )
}

expected_treatments <- c("Sham", "MCAO1hr", "MCAO3hr")
expected_sexes <- c("Female", "Male")
expected_ages <- c("Young", "Old")
if (
  !setequal(unique(manifest$treatment), expected_treatments) ||
    !setequal(unique(manifest$sex), expected_sexes) ||
    !setequal(unique(manifest$age_group), expected_ages)
) {
  stop(
    "Manifest treatment, sex, or age labels differ from expected values.",
    call. = FALSE
  )
}

manifest$treatment <- factor(
  manifest$treatment,
  levels = expected_treatments
)
manifest$sex <- factor(manifest$sex, levels = expected_sexes)
manifest$age_group <- factor(
  manifest$age_group,
  levels = expected_ages
)

covariate_design <- model.matrix(
  ~ treatment + sex + age_group,
  data = manifest
)
if (qr(covariate_design)$rank != ncol(covariate_design)) {
  stop(
    "The treatment + sex + age covariate design is not full-rank.",
    call. = FALSE
  )
}

if (!"ProbeSetName" %in% names(mirna_expression)) {
  stop(
    "miRNA expression matrix lacks ProbeSetName.",
    call. = FALSE
  )
}
if (!"transcript_cluster_id" %in% names(mrna_expression)) {
  stop(
    "mRNA expression matrix lacks transcript_cluster_id.",
    call. = FALSE
  )
}
if (
  anyDuplicated(as.character(mirna_expression$ProbeSetName)) ||
    anyDuplicated(as.character(mrna_expression$transcript_cluster_id))
) {
  stop(
    "An expression matrix contains duplicated feature identifiers.",
    call. = FALSE
  )
}
if (
  !all(manifest$mirna_file_name %in% names(mirna_expression)) ||
    !all(manifest$mrna_file_name %in% names(mrna_expression))
) {
  stop(
    "Not every paired-manifest sample occurs in both expression matrices.",
    call. = FALSE
  )
}

mirna_row <- match(
  as.character(candidates$mirna_probeset_id),
  as.character(mirna_expression$ProbeSetName)
)
mrna_row <- match(
  as.character(candidates$transcript_cluster_id),
  as.character(mrna_expression$transcript_cluster_id)
)
if (anyNA(mirna_row) || anyNA(mrna_row)) {
  stop(
    "Not every candidate feature was found in its expression matrix.",
    call. = FALSE
  )
}

association_rows <- vector("list", nrow(candidates))

for (i in seq_len(nrow(candidates))) {
  mirna_values <- as.numeric(unlist(
    mirna_expression[
      mirna_row[i],
      manifest$mirna_file_name,
      drop = FALSE
    ],
    use.names = FALSE
  ))
  mrna_values <- as.numeric(unlist(
    mrna_expression[
      mrna_row[i],
      manifest$mrna_file_name,
      drop = FALSE
    ],
    use.names = FALSE
  ))

  if (
    !all(is.finite(mirna_values)) ||
      !all(is.finite(mrna_values)) ||
      stats::sd(mirna_values) == 0 ||
      stats::sd(mrna_values) == 0
  ) {
    stop(
      "Pair ", pair_key[i],
      " has non-finite or constant expression values.",
      call. = FALSE
    )
  }

  analysis_data <- data.frame(
    mrna_expression = mrna_values,
    mirna_expression = mirna_values,
    treatment = manifest$treatment,
    sex = manifest$sex,
    age_group = manifest$age_group
  )

  fit <- stats::lm(
    mrna_expression ~ mirna_expression + treatment + sex + age_group,
    data = analysis_data
  )
  if (fit$rank != length(stats::coef(fit)) || anyNA(stats::coef(fit))) {
    stop(
      "Association model is rank-deficient for pair ", pair_key[i], ".",
      call. = FALSE
    )
  }

  coefficient_table <- summary(fit)$coefficients
  mirna_coefficient <- coefficient_table["mirna_expression", ]
  residual_df <- stats::df.residual(fit)
  critical_t <- stats::qt(0.975, df = residual_df)
  partial_r <- sign(mirna_coefficient["t value"]) * sqrt(
    mirna_coefficient["t value"]^2 /
      (mirna_coefficient["t value"]^2 + residual_df)
  )

  pearson_test <- stats::cor.test(
    mirna_values,
    mrna_values,
    method = "pearson"
  )

  association_rows[[i]] <- data.frame(
    n_paired_animals = stats::nobs(fit),
    unadjusted_pearson_r = unname(pearson_test$estimate),
    unadjusted_P.Value = pearson_test$p.value,
    adjusted_beta = unname(mirna_coefficient["Estimate"]),
    adjusted_SE = unname(mirna_coefficient["Std. Error"]),
    adjusted_t = unname(mirna_coefficient["t value"]),
    adjusted_residual_df = residual_df,
    adjusted_partial_r = unname(partial_r),
    adjusted_ci_lower = unname(
      mirna_coefficient["Estimate"] -
        critical_t * mirna_coefficient["Std. Error"]
    ),
    adjusted_ci_upper = unname(
      mirna_coefficient["Estimate"] +
        critical_t * mirna_coefficient["Std. Error"]
    ),
    adjusted_P.Value = unname(mirna_coefficient["Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
}

association_statistics <- do.call(rbind, association_rows)
association_statistics$adjusted_adj.P.Val <- stats::p.adjust(
  association_statistics$adjusted_P.Value,
  method = "BH"
)
association_statistics$association_direction <- ifelse(
  association_statistics$adjusted_beta < 0,
  "negative",
  "positive"
)
association_statistics$negative_adjusted_association <-
  association_statistics$adjusted_beta < 0
association_statistics$adjusted_fdr_lt_0_05 <-
  association_statistics$adjusted_adj.P.Val < 0.05
association_statistics$adjusted_fdr_lt_0_10 <-
  association_statistics$adjusted_adj.P.Val < 0.10

metadata_columns <- c(
  "mirna_id", "mirna_probeset_id", "transcript_cluster_id",
  "gene_name", "gene_entrez_id", "gene_symbol", "contrast",
  "mirna_logFC", "mirna_P.Value", "mirna_adj.P.Val",
  "mrna_logFC", "mrna_P.Value", "mrna_adj.P.Val",
  "inverse_logFC_direction",
  "mirtarbase_supported", "targetscanmouse_supported",
  "both_sources_supported", "evidence_sources"
)
results <- cbind(
  candidates[, metadata_columns, drop = FALSE],
  association_statistics
)
results$association_rank <- rank(
  results$adjusted_P.Value,
  ties.method = "min"
)

output_order <- order(
  results$adjusted_P.Value,
  results$mirna_id,
  results$gene_symbol
)
results <- results[output_order, , drop = FALSE]
row.names(results) <- NULL

if (
  nrow(results) != nrow(candidates) ||
    anyNA(results$adjusted_P.Value) ||
    anyNA(results$adjusted_adj.P.Val) ||
    any(!is.finite(results$adjusted_beta)) ||
    !all(results$n_paired_animals == nrow(manifest))
) {
  stop(
    "Association output failed final completeness checks.",
    call. = FALSE
  )
}

group_counts <- table(manifest$treatment)
top_result <- results[1L, , drop = FALSE]
summary_lines <- c(
  "Adjusted paired-animal miRNA--mRNA association summary",
  "",
  paste0("Candidate input: ", candidate_file),
  paste0("Paired manifest: ", manifest_file),
  paste0("miRNA expression input: ", mirna_expression_file),
  paste0("mRNA expression input: ", mrna_expression_file),
  paste0("Association result: ", result_file),
  "",
  paste0("Paired animals: ", nrow(manifest)),
  paste0(
    "Treatment counts: ",
    paste(names(group_counts), group_counts, sep = "=", collapse = ", ")
  ),
  paste0("Measured candidate pairs tested: ", nrow(results)),
  "Model: mRNA expression ~ miRNA expression + treatment + sex + age",
  "Association term tested: miRNA expression coefficient",
  "Multiple-testing correction: Benjamini-Hochberg across tested pairs",
  "",
  paste0(
    "Negative adjusted associations: ",
    sum(results$negative_adjusted_association)
  ),
  paste0(
    "Positive adjusted associations: ",
    sum(!results$negative_adjusted_association)
  ),
  paste0(
    "Adjusted association raw p-value < 0.05: ",
    sum(results$adjusted_P.Value < 0.05)
  ),
  paste0(
    "Adjusted association FDR < 0.10: ",
    sum(results$adjusted_fdr_lt_0_10)
  ),
  paste0(
    "Adjusted association FDR < 0.05: ",
    sum(results$adjusted_fdr_lt_0_05)
  ),
  "",
  paste0(
    "Smallest adjusted association p-value pair: ",
    top_result$mirna_id, " -- ", top_result$gene_symbol
  ),
  paste0(
    "Smallest adjusted association p-value: ",
    format(top_result$adjusted_P.Value, digits = 8)
  ),
  paste0(
    "Corresponding pair-level FDR: ",
    format(top_result$adjusted_adj.P.Val, digits = 8)
  ),
  paste0(
    "Adjusted association direction: ",
    top_result$association_direction
  ),
  "",
  paste(
    "Interpretation boundary: an adjusted expression association is",
    "supportive evidence only and does not establish direct regulation,",
    "causality, or mediation."
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, result_file, row.names = FALSE, na = "")
writeLines(summary_lines, summary_file)

cat("Paired animals: ", nrow(manifest), "\n", sep = "")
cat("Measured candidate pairs tested: ", nrow(results), "\n", sep = "")
cat(
  "Adjusted association raw p-value < 0.05: ",
  sum(results$adjusted_P.Value < 0.05),
  "\n",
  sep = ""
)
cat(
  "Adjusted association FDR < 0.10: ",
  sum(results$adjusted_fdr_lt_0_10),
  "\n",
  sep = ""
)
cat(
  "Adjusted association FDR < 0.05: ",
  sum(results$adjusted_fdr_lt_0_05),
  "\n",
  sep = ""
)
cat("Wrote: ", result_file, "\n", sep = "")
cat("Wrote: ", summary_file, "\n", sep = "")
