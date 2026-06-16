#!/usr/bin/env Rscript

# Examine how treatment, sex, age, and sample selection affect the apparent
# expression relationships among the 56 DE-supported miRNA--mRNA pairs.
#
# Five pre-specified models are compared:
#
#   1. mRNA ~ miRNA
#   2. mRNA ~ miRNA + treatment
#   3. mRNA ~ miRNA + sex + age
#   4. mRNA ~ miRNA + treatment + sex + age
#   5. Model 4 restricted to MCAO1hr and MCAO3hr animals
#
# Benjamini--Hochberg correction is applied separately across the 56 pairs
# within each model. The all-animal fully adjusted results are required to
# reproduce the primary results from Script 08 exactly.
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
#   results/multiomics/paired_association/
#     mirna_mrna_adjusted_association_results.csv
#
# Outputs:
#   results/multiomics/paired_association/covariate_effects/
#     mirna_mrna_association_model_comparison_long.csv
#     mirna_mrna_association_model_comparison_wide.csv
#     mirna_mrna_association_model_comparison_summary.txt
#     figures/
#       top_10_unadjusted_association_plots.png
#       inverse_de_pair_association_plots.png
#       unadjusted_vs_fully_adjusted_pvalues.png
#
# Usage:
#   Rscript \
#     scripts/multiomics/09_examine_covariate_effects_on_associations.R

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
primary_association_file <- file.path(
  "results", "multiomics", "paired_association",
  "mirna_mrna_adjusted_association_results.csv"
)
output_dir <- file.path(
  "results", "multiomics", "paired_association", "covariate_effects"
)
figure_dir <- file.path(output_dir, "figures")
long_output_file <- file.path(
  output_dir, "mirna_mrna_association_model_comparison_long.csv"
)
wide_output_file <- file.path(
  output_dir, "mirna_mrna_association_model_comparison_wide.csv"
)
summary_output_file <- file.path(
  output_dir, "mirna_mrna_association_model_comparison_summary.txt"
)
top_plot_file <- file.path(
  figure_dir, "top_10_unadjusted_association_plots.png"
)
inverse_plot_file <- file.path(
  figure_dir, "inverse_de_pair_association_plots.png"
)
pvalue_plot_file <- file.path(
  figure_dir, "unadjusted_vs_fully_adjusted_pvalues.png"
)

required_inputs <- c(
  candidate_file,
  manifest_file,
  mirna_expression_file,
  mrna_expression_file,
  primary_association_file
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
primary_associations <- read_table(primary_association_file)

required_candidate_columns <- c(
  "mirna_id", "mirna_probeset_id", "transcript_cluster_id",
  "gene_name", "gene_entrez_id", "gene_symbol", "contrast",
  "mirna_logFC", "mirna_P.Value", "mirna_adj.P.Val",
  "mrna_logFC", "mrna_P.Value", "mrna_adj.P.Val",
  "inverse_logFC_direction",
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

required_primary_columns <- c(
  "mirna_probeset_id", "transcript_cluster_id",
  "adjusted_beta", "adjusted_P.Value", "adjusted_adj.P.Val"
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

if (nrow(candidates) != 56L || nrow(manifest) != 43L) {
  stop(
    "Expected 56 candidate pairs and 43 paired animals; found ",
    nrow(candidates), " and ", nrow(manifest), ".",
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
primary_pair_key <- paste(
  primary_associations$mirna_probeset_id,
  primary_associations$transcript_cluster_id,
  sep = "||"
)
if (
  anyDuplicated(pair_key) ||
    anyDuplicated(primary_pair_key) ||
    !setequal(pair_key, primary_pair_key)
) {
  stop(
    "Candidate and primary-association measured-pair keys do not agree.",
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

if (!"ProbeSetName" %in% names(mirna_expression)) {
  stop("miRNA expression matrix lacks ProbeSetName.", call. = FALSE)
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
    "Not every paired animal occurs in both expression matrices.",
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

mirna_value_matrix <- data.matrix(
  mirna_expression[
    mirna_row,
    manifest$mirna_file_name,
    drop = FALSE
  ]
)
mrna_value_matrix <- data.matrix(
  mrna_expression[
    mrna_row,
    manifest$mrna_file_name,
    drop = FALSE
  ]
)
if (
  any(!is.finite(mirna_value_matrix)) ||
    any(!is.finite(mrna_value_matrix))
) {
  stop(
    "Candidate expression matrices contain non-finite values.",
    call. = FALSE
  )
}

mcao_sample_index <- which(
  manifest$treatment %in% c("MCAO1hr", "MCAO3hr")
)
if (
  length(mcao_sample_index) != 27L ||
    sum(manifest$treatment[mcao_sample_index] == "MCAO1hr") != 16L ||
    sum(manifest$treatment[mcao_sample_index] == "MCAO3hr") != 11L
) {
  stop(
    "The MCAO-focused subset does not contain the expected 27 animals.",
    call. = FALSE
  )
}

model_specs <- list(
  list(
    model_id = "unadjusted_all_43",
    model_label = "Unadjusted, all 43 animals",
    formula = mrna_expression ~ mirna_expression,
    sample_index = seq_len(nrow(manifest))
  ),
  list(
    model_id = "treatment_adjusted_all_43",
    model_label = "Treatment-adjusted, all 43 animals",
    formula = mrna_expression ~ mirna_expression + treatment,
    sample_index = seq_len(nrow(manifest))
  ),
  list(
    model_id = "sex_age_adjusted_all_43",
    model_label = "Sex+age-adjusted, all 43 animals",
    formula = mrna_expression ~ mirna_expression + sex + age_group,
    sample_index = seq_len(nrow(manifest))
  ),
  list(
    model_id = "fully_adjusted_all_43",
    model_label = "Treatment+sex+age-adjusted, all 43 animals",
    formula = (
      mrna_expression ~ mirna_expression + treatment + sex + age_group
    ),
    sample_index = seq_len(nrow(manifest))
  ),
  list(
    model_id = "fully_adjusted_mcao_27",
    model_label = "Treatment+sex+age-adjusted, MCAO 27 animals",
    formula = (
      mrna_expression ~ mirna_expression + treatment + sex + age_group
    ),
    sample_index = mcao_sample_index
  )
)
model_ids <- vapply(model_specs, `[[`, character(1), "model_id")

fit_pair_model <- function(
    mirna_values,
    mrna_values,
    sample_metadata,
    model_formula
) {
  analysis_data <- data.frame(
    mrna_expression = mrna_values,
    mirna_expression = mirna_values,
    treatment = droplevels(sample_metadata$treatment),
    sex = droplevels(sample_metadata$sex),
    age_group = droplevels(sample_metadata$age_group)
  )

  design <- model.matrix(model_formula, data = analysis_data)
  if (qr(design)$rank != ncol(design)) {
    stop("A requested association design is not full-rank.", call. = FALSE)
  }

  fit <- stats::lm(model_formula, data = analysis_data)
  coefficient_table <- summary(fit)$coefficients
  if (
    !"mirna_expression" %in% row.names(coefficient_table) ||
      anyNA(stats::coef(fit))
  ) {
    stop(
      "A requested association model lacks an estimable miRNA coefficient.",
      call. = FALSE
    )
  }

  mirna_coefficient <- coefficient_table["mirna_expression", ]
  residual_df <- stats::df.residual(fit)
  critical_t <- stats::qt(0.975, df = residual_df)
  partial_r <- sign(mirna_coefficient["t value"]) * sqrt(
    mirna_coefficient["t value"]^2 /
      (mirna_coefficient["t value"]^2 + residual_df)
  )

  data.frame(
    n_animals = stats::nobs(fit),
    beta = unname(mirna_coefficient["Estimate"]),
    standard_error = unname(mirna_coefficient["Std. Error"]),
    t_value = unname(mirna_coefficient["t value"]),
    residual_df = residual_df,
    partial_r = unname(partial_r),
    ci_lower = unname(
      mirna_coefficient["Estimate"] -
        critical_t * mirna_coefficient["Std. Error"]
    ),
    ci_upper = unname(
      mirna_coefficient["Estimate"] +
        critical_t * mirna_coefficient["Std. Error"]
    ),
    P.Value = unname(mirna_coefficient["Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
}

metadata_columns <- c(
  "mirna_id", "mirna_probeset_id", "transcript_cluster_id",
  "gene_name", "gene_entrez_id", "gene_symbol", "contrast",
  "mirna_logFC", "mirna_P.Value", "mirna_adj.P.Val",
  "mrna_logFC", "mrna_P.Value", "mrna_adj.P.Val",
  "inverse_logFC_direction",
  "mirtarbase_supported", "targetscanmouse_supported",
  "both_sources_supported", "evidence_sources"
)

model_results <- vector("list", length(model_specs))
for (model_number in seq_along(model_specs)) {
  model_spec <- model_specs[[model_number]]
  sample_index <- model_spec$sample_index
  sample_metadata <- manifest[sample_index, , drop = FALSE]

  pair_results <- vector("list", nrow(candidates))
  for (pair_number in seq_len(nrow(candidates))) {
    pair_results[[pair_number]] <- fit_pair_model(
      mirna_values = mirna_value_matrix[pair_number, sample_index],
      mrna_values = mrna_value_matrix[pair_number, sample_index],
      sample_metadata = sample_metadata,
      model_formula = model_spec$formula
    )
  }

  model_statistics <- do.call(rbind, pair_results)
  model_statistics$adj.P.Val <- stats::p.adjust(
    model_statistics$P.Value,
    method = "BH"
  )
  model_statistics$association_direction <- ifelse(
    model_statistics$beta < 0,
    "negative",
    "positive"
  )
  model_statistics$raw_p_lt_0_05 <- model_statistics$P.Value < 0.05
  model_statistics$fdr_lt_0_10 <- model_statistics$adj.P.Val < 0.10
  model_statistics$fdr_lt_0_05 <- model_statistics$adj.P.Val < 0.05
  model_statistics$association_rank <- rank(
    model_statistics$P.Value,
    ties.method = "min"
  )

  model_results[[model_number]] <- cbind(
    candidates[, metadata_columns, drop = FALSE],
    model_id = model_spec$model_id,
    model_label = model_spec$model_label,
    model_statistics
  )
}

long_results <- do.call(rbind, model_results)
long_results$model_id <- factor(
  long_results$model_id,
  levels = model_ids
)
long_results <- long_results[
  order(
    long_results$model_id,
    long_results$P.Value,
    long_results$mirna_id,
    long_results$gene_symbol
  ),
  ,
  drop = FALSE
]
long_results$model_id <- as.character(long_results$model_id)
row.names(long_results) <- NULL

wide_results <- candidates[, metadata_columns, drop = FALSE]
for (model_id in model_ids) {
  model_data <- long_results[
    long_results$model_id == model_id,
    ,
    drop = FALSE
  ]
  model_key <- paste(
    model_data$mirna_probeset_id,
    model_data$transcript_cluster_id,
    sep = "||"
  )
  model_match <- match(pair_key, model_key)
  if (anyNA(model_match)) {
    stop(
      "A model output is missing a candidate pair.",
      call. = FALSE
    )
  }

  for (statistic in c(
    "n_animals", "beta", "partial_r", "P.Value", "adj.P.Val",
    "association_direction", "raw_p_lt_0_05",
    "fdr_lt_0_10", "fdr_lt_0_05"
  )) {
    wide_results[[paste(model_id, statistic, sep = "_")]] <-
      model_data[[statistic]][model_match]
  }
}

wide_results$lost_raw_significance_after_treatment_adjustment <-
  wide_results$unadjusted_all_43_raw_p_lt_0_05 &
    !wide_results$treatment_adjusted_all_43_raw_p_lt_0_05
wide_results$lost_raw_significance_after_sex_age_adjustment <-
  wide_results$unadjusted_all_43_raw_p_lt_0_05 &
    !wide_results$sex_age_adjusted_all_43_raw_p_lt_0_05
wide_results$lost_raw_significance_after_full_adjustment <-
  wide_results$unadjusted_all_43_raw_p_lt_0_05 &
    !wide_results$fully_adjusted_all_43_raw_p_lt_0_05

primary_match <- match(pair_key, primary_pair_key)
full_match <- match(
  pair_key,
  paste(
    wide_results$mirna_probeset_id,
    wide_results$transcript_cluster_id,
    sep = "||"
  )
)
if (
  anyNA(primary_match) ||
    anyNA(full_match) ||
    !isTRUE(all.equal(
      wide_results$fully_adjusted_all_43_beta[full_match],
      primary_associations$adjusted_beta[primary_match],
      tolerance = 1e-12
    )) ||
    !isTRUE(all.equal(
      wide_results$fully_adjusted_all_43_P.Value[full_match],
      primary_associations$adjusted_P.Value[primary_match],
      tolerance = 1e-12
    )) ||
    !isTRUE(all.equal(
      wide_results$fully_adjusted_all_43_adj.P.Val[full_match],
      primary_associations$adjusted_adj.P.Val[primary_match],
      tolerance = 1e-12
    ))
) {
  stop(
    "The fully adjusted all-animal model does not reproduce Script 08.",
    call. = FALSE
  )
}

if (
  nrow(long_results) != nrow(candidates) * length(model_specs) ||
    anyNA(long_results$P.Value) ||
    anyNA(long_results$adj.P.Val) ||
    any(!is.finite(long_results$beta)) ||
    any(!is.finite(long_results$partial_r))
) {
  stop(
    "Model-comparison output failed completeness checks.",
    call. = FALSE
  )
}

model_summary <- do.call(
  rbind,
  lapply(model_ids, function(model_id) {
    model_data <- long_results[
      long_results$model_id == model_id,
      ,
      drop = FALSE
    ]
    top_row <- model_data[which.min(model_data$P.Value), , drop = FALSE]
    data.frame(
      model_id = model_id,
      model_label = unique(model_data$model_label),
      n_animals = unique(model_data$n_animals),
      tested_pairs = nrow(model_data),
      negative_associations = sum(
        model_data$association_direction == "negative"
      ),
      positive_associations = sum(
        model_data$association_direction == "positive"
      ),
      raw_p_lt_0_05 = sum(model_data$raw_p_lt_0_05),
      fdr_lt_0_10 = sum(model_data$fdr_lt_0_10),
      fdr_lt_0_05 = sum(model_data$fdr_lt_0_05),
      smallest_p_pair = paste(
        top_row$mirna_id,
        top_row$gene_symbol,
        sep = " -- "
      ),
      smallest_p_value = top_row$P.Value,
      smallest_p_fdr = top_row$adj.P.Val,
      stringsAsFactors = FALSE
    )
  })
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(long_results, long_output_file, row.names = FALSE, na = "")
write.csv(wide_results, wide_output_file, row.names = FALSE, na = "")

treatment_colors <- c(
  Sham = "#6B7280",
  MCAO1hr = "#D97706",
  MCAO3hr = "#2563EB"
)

plot_pair_panels <- function(pair_indices, output_file, plot_title) {
  n_panels <- length(pair_indices)
  n_columns <- if (n_panels <= 10L) 5L else 4L
  n_rows <- ceiling(n_panels / n_columns)

  grDevices::png(
    filename = output_file,
    width = 850L * n_columns,
    height = 700L * n_rows,
    res = 150
  )
  old_par <- graphics::par(
    mfrow = c(n_rows, n_columns),
    mar = c(5.4, 4.4, 4.0, 1.0),
    oma = c(0, 0, 2.5, 0)
  )
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  })

  for (panel_number in seq_along(pair_indices)) {
    pair_number <- pair_indices[panel_number]
    x <- mirna_value_matrix[pair_number, ]
    y <- mrna_value_matrix[pair_number, ]
    group <- manifest$treatment

    graphics::plot(
      x,
      y,
      pch = 16,
      col = grDevices::adjustcolor(
        treatment_colors[as.character(group)],
        alpha.f = 0.78
      ),
      xlab = "miRNA expression",
      ylab = "mRNA expression",
      main = paste0(
        candidates$mirna_id[pair_number],
        " \u2192 ",
        candidates$gene_symbol[pair_number]
      ),
      sub = paste0(
        "Unadjusted p=",
        format(
          wide_results$unadjusted_all_43_P.Value[pair_number],
          digits = 3
        ),
        "; full-adjusted p=",
        format(
          wide_results$fully_adjusted_all_43_P.Value[pair_number],
          digits = 3
        )
      )
    )

    graphics::abline(
      stats::lm(y ~ x),
      col = "#111827",
      lwd = 2,
      lty = 2
    )
    for (treatment_level in levels(group)) {
      use_group <- group == treatment_level
      if (
        sum(use_group) >= 3L &&
          stats::sd(x[use_group]) > 0
      ) {
        graphics::abline(
          stats::lm(y[use_group] ~ x[use_group]),
          col = treatment_colors[treatment_level],
          lwd = 2
        )
      }
    }

    if (panel_number == 1L) {
      graphics::legend(
        "topleft",
        legend = c(names(treatment_colors), "All animals"),
        col = c(treatment_colors, "#111827"),
        pch = c(rep(16, length(treatment_colors)), NA),
        lty = c(rep(1, length(treatment_colors)), 2),
        lwd = c(rep(2, length(treatment_colors)), 2),
        cex = 0.72,
        bg = grDevices::adjustcolor("white", alpha.f = 0.85)
      )
    }
  }

  unused_panels <- n_rows * n_columns - n_panels
  if (unused_panels > 0L) {
    for (unused in seq_len(unused_panels)) {
      graphics::plot.new()
    }
  }
  graphics::mtext(
    plot_title,
    outer = TRUE,
    side = 3,
    line = 0.5,
    font = 2,
    cex = 1.2
  )
}

top_pair_indices <- order(
  wide_results$unadjusted_all_43_P.Value
)[seq_len(10L)]
inverse_pair_indices <- which(candidates$inverse_logFC_direction)
if (length(inverse_pair_indices) != 13L) {
  stop(
    "Expected 13 inverse-DE candidate pairs; found ",
    length(inverse_pair_indices), ".",
    call. = FALSE
  )
}

plot_pair_panels(
  pair_indices = top_pair_indices,
  output_file = top_plot_file,
  plot_title = paste(
    "Ten smallest unadjusted association p-values;",
    "points and fitted lines are colored by treatment"
  )
)
plot_pair_panels(
  pair_indices = inverse_pair_indices,
  output_file = inverse_plot_file,
  plot_title = paste(
    "Thirteen inverse-DE candidate pairs;",
    "points and fitted lines are colored by treatment"
  )
)

grDevices::png(
  filename = pvalue_plot_file,
  width = 1800,
  height = 1500,
  res = 180
)
old_par <- graphics::par(mar = c(5.2, 5.3, 2.0, 1.2))
unadjusted_log_p <- -log10(
  wide_results$unadjusted_all_43_P.Value
)
fully_adjusted_log_p <- -log10(
  wide_results$fully_adjusted_all_43_P.Value
)
point_colors <- ifelse(
  candidates$inverse_logFC_direction,
  "#B91C1C",
  "#2563EB"
)
graphics::plot(
  unadjusted_log_p,
  fully_adjusted_log_p,
  pch = 16,
  col = grDevices::adjustcolor(point_colors, alpha.f = 0.75),
  xlim = c(
    max(0, min(unadjusted_log_p) - 0.05),
    max(unadjusted_log_p) + 0.35
  ),
  ylim = c(
    max(0, min(fully_adjusted_log_p) - 0.05),
    max(fully_adjusted_log_p) + 0.15
  ),
  xlab = expression(-log[10]("unadjusted p-value")),
  ylab = expression(-log[10]("fully adjusted p-value")),
  main = "Association evidence before and after covariate adjustment"
)
graphics::abline(
  a = 0,
  b = 1,
  lty = 2,
  col = "#4B5563"
)
graphics::abline(
  v = -log10(0.05),
  h = -log10(0.05),
  lty = 3,
  col = "#111827"
)
label_indices <- order(
  wide_results$unadjusted_all_43_P.Value
)[seq_len(5L)]
graphics::text(
  unadjusted_log_p[label_indices],
  fully_adjusted_log_p[label_indices],
  labels = paste(
    candidates$mirna_id[label_indices],
    candidates$gene_symbol[label_indices],
    sep = "\n"
  ),
  pos = 3,
  cex = 0.72
)
graphics::legend(
  "topleft",
  legend = c("Same DE direction", "Inverse DE direction"),
  col = c("#2563EB", "#B91C1C"),
  pch = 16,
  bty = "n"
)
graphics::par(old_par)
grDevices::dev.off()

summary_lines <- c(
  "Covariate effects on paired miRNA--mRNA associations",
  "",
  paste0("Candidate input: ", candidate_file),
  paste0("Primary association input: ", primary_association_file),
  paste0("All paired animals: ", nrow(manifest)),
  paste0("MCAO1hr/MCAO3hr focused animals: ", length(mcao_sample_index)),
  paste0("Candidate pairs tested per model: ", nrow(candidates)),
  "Multiple-testing correction: Benjamini-Hochberg within each model",
  "",
  "Model-level results:"
)
for (i in seq_len(nrow(model_summary))) {
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ", model_summary$model_label[i],
      ": raw p < 0.05 = ", model_summary$raw_p_lt_0_05[i],
      ", FDR < 0.10 = ", model_summary$fdr_lt_0_10[i],
      ", FDR < 0.05 = ", model_summary$fdr_lt_0_05[i],
      "; smallest-p pair = ", model_summary$smallest_p_pair[i],
      ", p = ",
      format(model_summary$smallest_p_value[i], digits = 7),
      ", FDR = ",
      format(model_summary$smallest_p_fdr[i], digits = 7)
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  paste0(
    "Unadjusted raw-p pairs lost after treatment adjustment: ",
    sum(
      wide_results$lost_raw_significance_after_treatment_adjustment
    )
  ),
  paste0(
    "Unadjusted raw-p pairs lost after sex+age adjustment: ",
    sum(
      wide_results$lost_raw_significance_after_sex_age_adjustment
    )
  ),
  paste0(
    "Unadjusted raw-p pairs lost after full adjustment: ",
    sum(
      wide_results$lost_raw_significance_after_full_adjustment
    )
  ),
  "",
  paste0("Wrote: ", long_output_file),
  paste0("Wrote: ", wide_output_file),
  paste0("Wrote: ", top_plot_file),
  paste0("Wrote: ", inverse_plot_file),
  paste0("Wrote: ", pvalue_plot_file),
  "",
  paste(
    "Interpretation boundary: model comparisons show whether apparent",
    "expression relationships persist after accounting for measured",
    "covariates. They do not establish direct regulation, causality,",
    "or mediation."
  )
)
writeLines(summary_lines, summary_output_file)

print(
  model_summary[
    ,
    c(
      "model_label", "n_animals", "raw_p_lt_0_05",
      "fdr_lt_0_10", "fdr_lt_0_05"
    )
  ],
  row.names = FALSE
)
cat(
  "Unadjusted raw-p pairs lost after treatment adjustment: ",
  sum(wide_results$lost_raw_significance_after_treatment_adjustment),
  "\n",
  sep = ""
)
cat(
  "Unadjusted raw-p pairs lost after sex+age adjustment: ",
  sum(wide_results$lost_raw_significance_after_sex_age_adjustment),
  "\n",
  sep = ""
)
cat(
  "Unadjusted raw-p pairs lost after full adjustment: ",
  sum(wide_results$lost_raw_significance_after_full_adjustment),
  "\n",
  sep = ""
)
cat("Wrote: ", long_output_file, "\n", sep = "")
cat("Wrote: ", wide_output_file, "\n", sep = "")
cat("Wrote: ", summary_output_file, "\n", sep = "")
cat("Wrote figures under: ", figure_dir, "\n", sep = "")
