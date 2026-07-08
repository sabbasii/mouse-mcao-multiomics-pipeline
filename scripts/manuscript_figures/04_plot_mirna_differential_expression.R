#!/usr/bin/env Rscript

# Create miRNA differential-expression panels for manuscript Figure 3.
#
# This script visualizes the saved covariate-adjusted miRNA limma results. It
# does not rerun differential expression. The displayed set is the union of
# reliably detected mature mouse miRNAs meeting BH FDR < 0.10 in at least one
# of the three prespecified treatment contrasts. For every selected miRNA, all
# three contrast estimates are shown.
#
# Approximate 95% confidence intervals are calculated from each saved
# moderated t statistic and log2 fold-change estimate. The resulting moderated
# standard error is combined with the residual degrees of freedom from the
# saved analysis design.
#
# Inputs:
#   results/mirna/differential_expression/rma_normalized_mirna/
#     treatment_sex_age_limma/
#       design_matrix_mirna.csv
#       dabg20_MCAO1hr_vs_Sham_mirna.csv
#       dabg20_MCAO3hr_vs_Sham_mirna.csv
#       dabg20_MCAO3hr_vs_MCAO1hr_mirna.csv
#
# Additional Figure 3B inputs:
#   results/mirna/expression/rma_normalized_mirna/annotation/
#     mouse_mature_mirna_expression_mirna.csv
#   results/mirna/differential_expression/rma_normalized_mirna/
#     treatment_sex_age_limma/analysis_samples_mirna.csv
#
# Outputs:
#   results/manuscript/figures/
#     components/figure3a_mirna_cross_contrast_effects.png
#     components/figure3b_mirna_6968_adjusted_expression.png
#     components/figure3c_mirna_8101_adjusted_expression.png
#     components/figure3d_mirna_5130_adjusted_expression.png
#     main/figure3_mirna_differential_expression.png
#     source_data/figure3a_mirna_cross_contrast_effects.csv
#     source_data/figure3b_mirna_6968_sample_expression.csv
#     source_data/figure3b_mirna_6968_adjusted_estimates.csv
#     source_data/figure3c_mirna_8101_sample_expression.csv
#     source_data/figure3c_mirna_8101_adjusted_estimates.csv
#     source_data/figure3d_mirna_5130_sample_expression.csv
#     source_data/figure3d_mirna_5130_adjusted_estimates.csv
#
# Usage:
#   Rscript scripts/manuscript_figures/04_plot_mirna_differential_expression.R

required_packages <- c("ggplot2", "limma")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

input_dir <- file.path(
  "results",
  "mirna",
  "differential_expression",
  "rma_normalized_mirna",
  "treatment_sex_age_limma"
)
design_file <- file.path(input_dir, "design_matrix_mirna.csv")
analysis_sample_file <- file.path(
  input_dir,
  "analysis_samples_mirna.csv"
)
expression_file <- file.path(
  "results",
  "mirna",
  "expression",
  "rma_normalized_mirna",
  "annotation",
  "mouse_mature_mirna_expression_mirna.csv"
)
figure_root <- file.path("results", "manuscript", "figures")
main_output_dir <- file.path(figure_root, "main")
component_output_dir <- file.path(figure_root, "components")
source_data_output_dir <- file.path(figure_root, "source_data")
figure_file <- file.path(
  component_output_dir,
  "figure3a_mirna_cross_contrast_effects.png"
)
source_data_file <- file.path(
  source_data_output_dir,
  "figure3a_mirna_cross_contrast_effects.csv"
)
expression_figure_file <- file.path(
  component_output_dir,
  "figure3b_mirna_6968_adjusted_expression.png"
)
sample_expression_file <- file.path(
  source_data_output_dir,
  "figure3b_mirna_6968_sample_expression.csv"
)
adjusted_estimate_file <- file.path(
  source_data_output_dir,
  "figure3b_mirna_6968_adjusted_estimates.csv"
)
mirna_8101_figure_file <- file.path(
  component_output_dir,
  "figure3c_mirna_8101_adjusted_expression.png"
)
mirna_8101_sample_file <- file.path(
  source_data_output_dir,
  "figure3c_mirna_8101_sample_expression.csv"
)
mirna_8101_adjusted_file <- file.path(
  source_data_output_dir,
  "figure3c_mirna_8101_adjusted_estimates.csv"
)
mirna_5130_figure_file <- file.path(
  component_output_dir,
  "figure3d_mirna_5130_adjusted_expression.png"
)
mirna_5130_sample_file <- file.path(
  source_data_output_dir,
  "figure3d_mirna_5130_sample_expression.csv"
)
mirna_5130_adjusted_file <- file.path(
  source_data_output_dir,
  "figure3d_mirna_5130_adjusted_estimates.csv"
)
main_figure_file <- file.path(
  main_output_dir,
  "figure3_mirna_differential_expression.png"
)

dir.create(
  main_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  component_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  source_data_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

contrast_order <- c(
  "MCAO1hr_vs_Sham",
  "MCAO3hr_vs_Sham",
  "MCAO3hr_vs_MCAO1hr"
)
contrast_labels <- c(
  MCAO1hr_vs_Sham = "MCAO1hr vs Sham",
  MCAO3hr_vs_Sham = "MCAO3hr vs Sham",
  MCAO3hr_vs_MCAO1hr = "MCAO3hr vs MCAO1hr"
)
input_files <- stats::setNames(
  file.path(
    input_dir,
    paste0("dabg20_", contrast_order, "_mirna.csv")
  ),
  contrast_order
)

required_inputs <- c(
  design_file,
  analysis_sample_file,
  expression_file,
  input_files
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required Figure 3A input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

design <- read.csv(
  design_file,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
design_matrix <- as.matrix(design)
storage.mode(design_matrix) <- "numeric"
if (
  nrow(design_matrix) != 43L ||
    ncol(design_matrix) != 5L ||
    any(!is.finite(design_matrix)) ||
    qr(design_matrix)$rank != ncol(design_matrix)
) {
  stop("Saved miRNA design matrix failed validation.", call. = FALSE)
}
residual_degrees_freedom <-
  nrow(design_matrix) - qr(design_matrix)$rank
critical_value <- stats::qt(
  0.975,
  df = residual_degrees_freedom
)

required_result_columns <- c(
  "ProbeSetName",
  "logFC",
  "t",
  "P.Value",
  "adj.P.Val",
  "Accession",
  "Transcript ID(Array Design)"
)
result_tables <- vector("list", length(contrast_order))
names(result_tables) <- contrast_order
reference_feature_ids <- NULL

for (contrast_name in contrast_order) {
  result <- read.csv(
    input_files[[contrast_name]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  missing_columns <- setdiff(required_result_columns, names(result))
  if (length(missing_columns) > 0L) {
    stop(
      basename(input_files[[contrast_name]]),
      " lacks required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (
    nrow(result) != 467L ||
      anyDuplicated(result$ProbeSetName) ||
      anyDuplicated(result[["Transcript ID(Array Design)"]]) ||
      any(!is.finite(result$logFC)) ||
      any(!is.finite(result$t)) ||
      any(!is.finite(result$P.Value)) ||
      any(!is.finite(result$adj.P.Val)) ||
      any(result$P.Value < 0 | result$P.Value > 1) ||
      any(result$adj.P.Val < 0 | result$adj.P.Val > 1)
  ) {
    stop(
      basename(input_files[[contrast_name]]),
      " failed saved-result validation.",
      call. = FALSE
    )
  }

  feature_ids <- sort(as.character(result$ProbeSetName))
  if (is.null(reference_feature_ids)) {
    reference_feature_ids <- feature_ids
  } else if (!identical(feature_ids, reference_feature_ids)) {
    stop(
      "The three reliably detected miRNA result sets do not contain ",
      "the same measured features.",
      call. = FALSE
    )
  }

  result$contrast <- contrast_name
  result_tables[[contrast_name]] <- result
}

all_results <- do.call(rbind, result_tables)
selected_feature_ids <- unique(
  all_results$ProbeSetName[all_results$adj.P.Val < 0.10]
)
if (
  length(selected_feature_ids) != 13L ||
    !all(
      c("mmu-miR-6968-5p", "mmu-miR-8101", "mmu-miR-5130") %in%
        all_results[["Transcript ID(Array Design)"]][
          all_results$ProbeSetName %in% selected_feature_ids
        ]
    )
) {
  stop(
    "The expected 13-miRNA FDR < 0.10 union was not recovered.",
    call. = FALSE
  )
}

plot_data <- all_results[
  all_results$ProbeSetName %in% selected_feature_ids,
  ,
  drop = FALSE
]
if (
  nrow(plot_data) != 39L ||
    any(plot_data$t == 0)
) {
  stop(
    "Selected cross-contrast miRNA results failed validation.",
    call. = FALSE
  )
}

plot_data$mirna_id <- sub(
  "^mmu-",
  "",
  plot_data[["Transcript ID(Array Design)"]]
)
plot_data$moderated_standard_error <- abs(
  plot_data$logFC / plot_data$t
)
plot_data$ci_lower <- plot_data$logFC -
  critical_value * plot_data$moderated_standard_error
plot_data$ci_upper <- plot_data$logFC +
  critical_value * plot_data$moderated_standard_error

plot_data$fdr_category <- ifelse(
  plot_data$adj.P.Val < 0.05,
  "FDR < 0.05",
  ifelse(
    plot_data$adj.P.Val < 0.10,
    "0.05 ≤ FDR < 0.10",
    "FDR ≥ 0.10"
  )
)
fdr_category_order <- c(
  "FDR < 0.05",
  "0.05 ≤ FDR < 0.10",
  "FDR ≥ 0.10"
)
fdr_colors <- c(
  "FDR < 0.05" = "#C8323E",
  "0.05 ≤ FDR < 0.10" = "#E69F00",
  "FDR ≥ 0.10" = "#A8B1BA"
)
plot_data$fdr_category <- factor(
  plot_data$fdr_category,
  levels = fdr_category_order
)

recovery_rows <- plot_data$contrast == "MCAO3hr_vs_MCAO1hr"
mirna_order <- plot_data$mirna_id[recovery_rows][
  order(
    plot_data$logFC[recovery_rows],
    decreasing = TRUE
  )
]
if (
  length(mirna_order) != 13L ||
    anyDuplicated(mirna_order)
) {
  stop("Could not derive a unique miRNA display order.", call. = FALSE)
}
plot_data$mirna_id <- factor(
  plot_data$mirna_id,
  levels = rev(mirna_order)
)
plot_data$contrast <- factor(
  plot_data$contrast,
  levels = contrast_order
)
plot_data$contrast_label <- factor(
  contrast_labels[as.character(plot_data$contrast)],
  levels = unname(contrast_labels[contrast_order])
)

source_output <- plot_data[
  order(
    plot_data$contrast,
    match(as.character(plot_data$mirna_id), rev(mirna_order))
  ),
  c(
    "ProbeSetName",
    "Accession",
    "Transcript ID(Array Design)",
    "mirna_id",
    "contrast",
    "logFC",
    "moderated_standard_error",
    "ci_lower",
    "ci_upper",
    "t",
    "P.Value",
    "adj.P.Val",
    "fdr_category"
  ),
  drop = FALSE
]
source_output$mirna_id <- as.character(source_output$mirna_id)
source_output$contrast <- as.character(source_output$contrast)
source_output$fdr_category <- as.character(
  source_output$fdr_category
)
write.csv(
  source_output,
  source_data_file,
  row.names = FALSE
)

coefficient_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = logFC,
    y = mirna_id,
    color = fdr_category
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "#455563",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = ci_lower,
      xend = ci_upper,
      yend = mirna_id
    ),
    linewidth = 1.15,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    size = 5.0,
    shape = 16
  ) +
  ggplot2::facet_grid(
    cols = ggplot2::vars(contrast_label)
  ) +
  ggplot2::scale_color_manual(
    values = fdr_colors,
    breaks = fdr_category_order,
    drop = FALSE,
    name = "Statistical evidence"
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(-0.5, 0, 0.5, 1.0),
    limits = c(-0.85, 1.15),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::labs(
    title = "(A) Recovery-associated miRNA effects across contrasts",
    subtitle = paste0(
      "Thirteen reliably detected miRNAs meeting FDR < 0.10 ",
      "in at least one comparison"
    ),
    x = expression(
      "Adjusted log"[2] * " fold change (95% CI)"
    ),
    y = NULL,
    caption = paste0(
      "Positive values indicate higher expression in the first group ",
      "named. Model: treatment + sex + age."
    )
  ) +
  ggplot2::theme_minimal(
    base_size = 18,
    base_family = "sans"
  ) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 28,
      color = "#17222C",
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = ggplot2::element_text(
      face = "bold",
      size = 18,
      color = "#50606D",
      margin = ggplot2::margin(b = 14)
    ),
    plot.caption = ggplot2::element_text(
      face = "bold",
      size = 15,
      color = "#43515D",
      hjust = 0,
      margin = ggplot2::margin(t = 12)
    ),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 17,
      color = "#26343E"
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 18,
      color = "#1E2B34",
      margin = ggplot2::margin(r = 8)
    ),
    axis.title.x = ggplot2::element_text(
      face = "bold",
      size = 19,
      color = "#1E2B34",
      margin = ggplot2::margin(t = 10)
    ),
    panel.grid.major.x = ggplot2::element_line(
      color = "#D9E0E5",
      linewidth = 0.45
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#EFF2F4",
      linewidth = 0.4
    ),
    panel.spacing.x = grid::unit(0.8, "cm"),
    strip.background = ggplot2::element_rect(
      fill = "#294C60",
      color = NA
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 18,
      color = "white",
      margin = ggplot2::margin(8, 8, 8, 8)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 16,
      color = "#26343E"
    ),
    legend.text = ggplot2::element_text(
      face = "bold",
      size = 15,
      color = "#34424D"
    ),
    legend.key.width = grid::unit(0.8, "cm"),
    plot.margin = ggplot2::margin(16, 18, 12, 16)
  )

grDevices::png(
  filename = figure_file,
  width = 14.0,
  height = 9.4,
  units = "in",
  res = 400,
  bg = "white"
)
print(coefficient_plot)
grDevices::dev.off()

analysis_samples <- read.csv(
  analysis_sample_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_sample_columns <- c(
  "file_name",
  "animal_id",
  "treatment",
  "sex",
  "age"
)
missing_sample_columns <- setdiff(
  required_sample_columns,
  names(analysis_samples)
)
if (
  length(missing_sample_columns) > 0L ||
    nrow(analysis_samples) != 43L ||
    anyDuplicated(analysis_samples$file_name) ||
    anyDuplicated(analysis_samples$animal_id)
) {
  stop(
    "Saved miRNA analysis samples failed validation.",
    call. = FALSE
  )
}

treatment_order <- c("Sham", "MCAO1hr", "MCAO3hr")
analysis_samples$treatment <- factor(
  analysis_samples$treatment,
  levels = treatment_order
)
analysis_samples$sex <- factor(
  analysis_samples$sex,
  levels = c("Female", "Male")
)
analysis_samples$age <- factor(
  analysis_samples$age,
  levels = c("Old", "Young")
)
if (
  anyNA(analysis_samples$treatment) ||
    anyNA(analysis_samples$sex) ||
    anyNA(analysis_samples$age)
) {
  stop(
    "Treatment, sex, or age is missing from the saved miRNA samples.",
    call. = FALSE
  )
}

rebuilt_design <- stats::model.matrix(
  ~ 0 + treatment + sex + age,
  data = analysis_samples
)
colnames(rebuilt_design) <- sub(
  "^treatment",
  "",
  colnames(rebuilt_design)
)
if (
  !identical(dim(design_matrix), dim(rebuilt_design)) ||
    !identical(colnames(design_matrix), colnames(rebuilt_design)) ||
    !isTRUE(
      all.equal(
        as.numeric(design_matrix),
        as.numeric(rebuilt_design),
        tolerance = 1e-12
      )
    )
) {
  stop(
    "The saved miRNA design does not reproduce treatment + sex + age.",
    call. = FALSE
  )
}

primary_statistics <- plot_data[
  as.character(plot_data$mirna_id) == "miR-6968-5p" &
    as.character(plot_data$contrast) == "MCAO3hr_vs_Sham",
  ,
  drop = FALSE
]
if (nrow(primary_statistics) != 1L) {
  stop(
    "Could not recover the primary miR-6968-5p result.",
    call. = FALSE
  )
}
target_probeset <- as.character(primary_statistics$ProbeSetName)

expression_table <- read.csv(
  expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (
  !"ProbeSetName" %in% names(expression_table) ||
    anyDuplicated(expression_table$ProbeSetName) ||
    !all(analysis_samples$file_name %in% names(expression_table))
) {
  stop(
    "The saved mature-miRNA expression matrix failed validation.",
    call. = FALSE
  )
}
expression_row <- match(
  target_probeset,
  as.character(expression_table$ProbeSetName)
)
if (is.na(expression_row)) {
  stop(
    "The miR-6968-5p probeset is absent from the expression matrix.",
    call. = FALSE
  )
}

selected_expression <- matrix(
  as.numeric(
    expression_table[
      expression_row,
      analysis_samples$file_name,
      drop = TRUE
    ]
  ),
  nrow = 1L,
  dimnames = list(
    "miR-6968-5p",
    analysis_samples$animal_id
  )
)
if (any(!is.finite(selected_expression))) {
  stop(
    "miR-6968-5p expression contains non-finite values.",
    call. = FALSE
  )
}

selected_fit <- limma::lmFit(
  selected_expression,
  design_matrix
)
coefficient <- selected_fit$coefficients[1L, ]
fitted_logfc <- coefficient[["MCAO3hr"]] -
  coefficient[["Sham"]]
if (
  !isTRUE(
    all.equal(
      unname(fitted_logfc),
      primary_statistics$logFC,
      tolerance = 1e-8
    )
  )
) {
  stop(
    "The selected-feature fit does not reproduce the saved ",
    "miR-6968-5p log2 fold change.",
    call. = FALSE
  )
}

sample_expression <- data.frame(
  probeset_id = target_probeset,
  mirna_id = "miR-6968-5p",
  animal_id = analysis_samples$animal_id,
  file_name = analysis_samples$file_name,
  treatment = analysis_samples$treatment,
  sex = analysis_samples$sex,
  age = analysis_samples$age,
  expression = as.numeric(selected_expression[1L, ]),
  stringsAsFactors = FALSE
)
sample_expression$treatment <- factor(
  sample_expression$treatment,
  levels = treatment_order
)

covariate_columns <- c("sexMale", "ageYoung")
if (!all(covariate_columns %in% colnames(design_matrix))) {
  stop(
    "The saved miRNA model lacks the expected sex and age covariates.",
    call. = FALSE
  )
}
covariate_means <- colMeans(
  design_matrix[, covariate_columns, drop = FALSE]
)

adjusted_rows <- vector("list", length(treatment_order))
for (index in seq_along(treatment_order)) {
  treatment_name <- treatment_order[[index]]
  estimate_vector <- stats::setNames(
    rep(0, ncol(design_matrix)),
    colnames(design_matrix)
  )
  estimate_vector[[treatment_name]] <- 1
  estimate_vector[covariate_columns] <- covariate_means

  estimate <- sum(coefficient * estimate_vector)
  unscaled_standard_error <- sqrt(
    drop(
      t(estimate_vector) %*%
        selected_fit$cov.coefficients %*%
        estimate_vector
    )
  )
  standard_error <- selected_fit$sigma[[1L]] *
    unscaled_standard_error
  degrees_freedom <- selected_fit$df.residual[[1L]]
  group_critical_value <- stats::qt(
    0.975,
    df = degrees_freedom
  )

  adjusted_rows[[index]] <- data.frame(
    probeset_id = target_probeset,
    mirna_id = "miR-6968-5p",
    treatment = treatment_name,
    adjusted_estimate = estimate,
    standard_error = standard_error,
    confidence_low = estimate -
      group_critical_value * standard_error,
    confidence_high = estimate +
      group_critical_value * standard_error,
    model = "expression ~ treatment + sex + age",
    stringsAsFactors = FALSE
  )
}
adjusted_estimates <- do.call(rbind, adjusted_rows)
adjusted_estimates$treatment <- factor(
  adjusted_estimates$treatment,
  levels = treatment_order
)

write.csv(
  sample_expression,
  sample_expression_file,
  row.names = FALSE
)
write.csv(
  adjusted_estimates,
  adjusted_estimate_file,
  row.names = FALSE
)

treatment_colors <- c(
  Sham = "#13A62A",
  MCAO1hr = "#1749D1",
  MCAO3hr = "#E32636"
)
treatment_counts <- table(sample_expression$treatment)
treatment_axis_labels <- paste0(
  treatment_order,
  "\n(n = ",
  as.integer(treatment_counts[treatment_order]),
  ")"
)

expression_range <- range(
  c(
    sample_expression$expression,
    adjusted_estimates$confidence_low,
    adjusted_estimates$confidence_high
  )
)
expression_span <- diff(expression_range)
if (!is.finite(expression_span) || expression_span <= 0) {
  stop("miR-6968-5p expression has an invalid range.", call. = FALSE)
}
bracket_y <- expression_range[[2L]] + 0.10 * expression_span
bracket_tip_y <- bracket_y - 0.025 * expression_span
bracket_label_y <- bracket_y + 0.045 * expression_span
plot_upper <- bracket_label_y + 0.08 * expression_span
plot_lower <- expression_range[[1L]] - 0.06 * expression_span

format_p_value <- function(value) {
  if (value < 0.001) {
    exponent <- floor(log10(value))
    mantissa <- value / (10^exponent)
    superscript_characters <- c(
      "-" = "⁻",
      "0" = "⁰",
      "1" = "¹",
      "2" = "²",
      "3" = "³",
      "4" = "⁴",
      "5" = "⁵",
      "6" = "⁶",
      "7" = "⁷",
      "8" = "⁸",
      "9" = "⁹"
    )
    exponent_characters <- strsplit(
      as.character(exponent),
      "",
      fixed = TRUE
    )[[1L]]
    superscript_exponent <- paste0(
      superscript_characters[exponent_characters],
      collapse = ""
    )
    paste0(
      sprintf("%.2f", mantissa),
      " × 10",
      superscript_exponent
    )
  } else {
    sprintf("%.3f", value)
  }
}

expression_plot <- ggplot2::ggplot(
  sample_expression,
  ggplot2::aes(
    x = treatment,
    y = expression,
    fill = treatment
  )
) +
  ggplot2::geom_violin(
    width = 0.72,
    scale = "width",
    trim = FALSE,
    alpha = 0.11,
    color = NA
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.12,
      height = 0,
      seed = 20260729
    ),
    shape = 21,
    size = 3.0,
    stroke = 0.38,
    color = "white",
    alpha = 0.82
  ) +
  ggplot2::geom_line(
    data = adjusted_estimates,
    ggplot2::aes(
      x = treatment,
      y = adjusted_estimate,
      group = 1
    ),
    inherit.aes = FALSE,
    linewidth = 0.85,
    color = "#34424C"
  ) +
  ggplot2::geom_errorbar(
    data = adjusted_estimates,
    ggplot2::aes(
      x = treatment,
      ymin = confidence_low,
      ymax = confidence_high
    ),
    inherit.aes = FALSE,
    width = 0.13,
    linewidth = 0.85,
    color = "#202B33"
  ) +
  ggplot2::geom_point(
    data = adjusted_estimates,
    ggplot2::aes(
      x = treatment,
      y = adjusted_estimate,
      fill = treatment
    ),
    inherit.aes = FALSE,
    shape = 23,
    size = 5.1,
    stroke = 0.95,
    color = "#202B33"
  ) +
  ggplot2::annotate(
    "segment",
    x = 1,
    xend = 3,
    y = bracket_y,
    yend = bracket_y,
    linewidth = 0.75,
    color = "#34424C"
  ) +
  ggplot2::annotate(
    "segment",
    x = c(1, 3),
    xend = c(1, 3),
    y = bracket_y,
    yend = bracket_tip_y,
    linewidth = 0.75,
    color = "#34424C"
  ) +
  ggplot2::annotate(
    "text",
    x = 2,
    y = bracket_label_y,
    label = paste0(
      "MCAO3hr vs Sham: FDR = ",
      sprintf("%.3f", primary_statistics$adj.P.Val)
    ),
    fontface = "bold",
    size = 5.2,
    color = "#283640"
  ) +
  ggplot2::scale_fill_manual(
    values = treatment_colors,
    guide = "none"
  ) +
  ggplot2::scale_x_discrete(
    labels = treatment_axis_labels
  ) +
  ggplot2::coord_cartesian(
    ylim = c(plot_lower, plot_upper),
    clip = "off"
  ) +
  ggplot2::labs(
    title = "(B) miR-6968-5p expression",
    subtitle = paste0(
      "MCAO3hr vs Sham\nAdjusted log2FC = ",
      sprintf("%+.2f", primary_statistics$logFC),
      "; p = ",
      format_p_value(primary_statistics$P.Value),
      "; FDR = ",
      sprintf("%.3f", primary_statistics$adj.P.Val)
    ),
    x = NULL,
    y = expression(log[2] * "-normalized expression"),
    caption = paste0(
      "Circles: individual animals. Diamonds/error bars: adjusted ",
      "means and 95% CIs."
    )
  ) +
  ggplot2::theme_minimal(
    base_size = 18,
    base_family = "sans"
  ) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 26,
      color = "#17222C",
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = ggplot2::element_text(
      face = "bold",
      size = 15.5,
      color = "#50606D",
      margin = ggplot2::margin(b = 12)
    ),
    plot.caption = ggplot2::element_text(
      face = "bold",
      size = 13,
      color = "#43515D",
      hjust = 0,
      margin = ggplot2::margin(t = 10)
    ),
    axis.title.y = ggplot2::element_text(
      face = "bold",
      size = 17,
      color = "#25323C",
      margin = ggplot2::margin(r = 10)
    ),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 16.5,
      color = "#25323C",
      lineheight = 0.95,
      margin = ggplot2::margin(t = 7)
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 16,
      color = "#34414D"
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E3E8EC",
      linewidth = 0.45
    ),
    panel.border = ggplot2::element_rect(
      color = "#5B6872",
      fill = NA,
      linewidth = 0.65
    ),
    plot.margin = ggplot2::margin(14, 16, 12, 14)
  )

grDevices::png(
  filename = expression_figure_file,
  width = 6.4,
  height = 6.7,
  units = "in",
  res = 400,
  bg = "white"
)
print(expression_plot)
grDevices::dev.off()

make_additional_expression_panel <- function(
  mirna_id,
  panel_label,
  highlighted_contrasts,
  output_figure_file,
  output_sample_file,
  output_adjusted_file,
  jitter_seed
) {
  selected_statistics <- do.call(
    rbind,
    lapply(
      highlighted_contrasts,
      function(contrast_name) {
        selected <- plot_data[
          as.character(plot_data$mirna_id) == mirna_id &
            as.character(plot_data$contrast) == contrast_name,
          ,
          drop = FALSE
        ]
        if (nrow(selected) != 1L) {
          stop(
            "Could not recover ",
            mirna_id,
            " for ",
            contrast_name,
            ".",
            call. = FALSE
          )
        }
        selected
      }
    )
  )
  target_probesets <- unique(
    as.character(selected_statistics$ProbeSetName)
  )
  if (length(target_probesets) != 1L) {
    stop(
      "Highlighted contrasts do not resolve to one probeset for ",
      mirna_id,
      ".",
      call. = FALSE
    )
  }
  target_probeset_local <- target_probesets[[1L]]
  expression_row_local <- match(
    target_probeset_local,
    as.character(expression_table$ProbeSetName)
  )
  if (is.na(expression_row_local)) {
    stop(
      mirna_id,
      " is absent from the expression matrix.",
      call. = FALSE
    )
  }

  selected_expression_local <- matrix(
    as.numeric(
      expression_table[
        expression_row_local,
        analysis_samples$file_name,
        drop = TRUE
      ]
    ),
    nrow = 1L,
    dimnames = list(
      mirna_id,
      analysis_samples$animal_id
    )
  )
  if (any(!is.finite(selected_expression_local))) {
    stop(
      mirna_id,
      " expression contains non-finite values.",
      call. = FALSE
    )
  }

  selected_fit_local <- limma::lmFit(
    selected_expression_local,
    design_matrix
  )
  coefficient_local <- selected_fit_local$coefficients[1L, ]
  fitted_contrast <- function(contrast_name) {
    if (contrast_name == "MCAO1hr_vs_Sham") {
      coefficient_local[["MCAO1hr"]] -
        coefficient_local[["Sham"]]
    } else if (contrast_name == "MCAO3hr_vs_Sham") {
      coefficient_local[["MCAO3hr"]] -
        coefficient_local[["Sham"]]
    } else if (contrast_name == "MCAO3hr_vs_MCAO1hr") {
      coefficient_local[["MCAO3hr"]] -
        coefficient_local[["MCAO1hr"]]
    } else {
      stop("Unsupported treatment contrast.", call. = FALSE)
    }
  }
  for (index in seq_len(nrow(selected_statistics))) {
    reproduced_logfc <- fitted_contrast(
      as.character(selected_statistics$contrast[[index]])
    )
    if (
      !isTRUE(
        all.equal(
          unname(reproduced_logfc),
          selected_statistics$logFC[[index]],
          tolerance = 1e-8
        )
      )
    ) {
      stop(
        "The selected-feature fit does not reproduce the saved ",
        mirna_id,
        " log2 fold change.",
        call. = FALSE
      )
    }
  }

  sample_expression_local <- data.frame(
    probeset_id = target_probeset_local,
    mirna_id = mirna_id,
    animal_id = analysis_samples$animal_id,
    file_name = analysis_samples$file_name,
    treatment = analysis_samples$treatment,
    sex = analysis_samples$sex,
    age = analysis_samples$age,
    expression = as.numeric(selected_expression_local[1L, ]),
    stringsAsFactors = FALSE
  )
  sample_expression_local$treatment <- factor(
    sample_expression_local$treatment,
    levels = treatment_order
  )

  adjusted_rows_local <- vector(
    "list",
    length(treatment_order)
  )
  for (index in seq_along(treatment_order)) {
    treatment_name <- treatment_order[[index]]
    estimate_vector <- stats::setNames(
      rep(0, ncol(design_matrix)),
      colnames(design_matrix)
    )
    estimate_vector[[treatment_name]] <- 1
    estimate_vector[covariate_columns] <- covariate_means

    estimate <- sum(coefficient_local * estimate_vector)
    unscaled_standard_error <- sqrt(
      drop(
        t(estimate_vector) %*%
          selected_fit_local$cov.coefficients %*%
          estimate_vector
      )
    )
    standard_error <- selected_fit_local$sigma[[1L]] *
      unscaled_standard_error
    degrees_freedom <- selected_fit_local$df.residual[[1L]]
    critical_value_local <- stats::qt(
      0.975,
      df = degrees_freedom
    )

    adjusted_rows_local[[index]] <- data.frame(
      probeset_id = target_probeset_local,
      mirna_id = mirna_id,
      treatment = treatment_name,
      adjusted_estimate = estimate,
      standard_error = standard_error,
      confidence_low = estimate -
        critical_value_local * standard_error,
      confidence_high = estimate +
        critical_value_local * standard_error,
      model = "expression ~ treatment + sex + age",
      stringsAsFactors = FALSE
    )
  }
  adjusted_estimates_local <- do.call(
    rbind,
    adjusted_rows_local
  )
  adjusted_estimates_local$treatment <- factor(
    adjusted_estimates_local$treatment,
    levels = treatment_order
  )

  write.csv(
    sample_expression_local,
    output_sample_file,
    row.names = FALSE
  )
  write.csv(
    adjusted_estimates_local,
    output_adjusted_file,
    row.names = FALSE
  )

  treatment_counts_local <- table(
    sample_expression_local$treatment
  )
  treatment_axis_labels_local <- paste0(
    treatment_order,
    "\n(n = ",
    as.integer(treatment_counts_local[treatment_order]),
    ")"
  )

  expression_range_local <- range(
    c(
      sample_expression_local$expression,
      adjusted_estimates_local$confidence_low,
      adjusted_estimates_local$confidence_high
    )
  )
  expression_span_local <- diff(expression_range_local)
  if (
    !is.finite(expression_span_local) ||
      expression_span_local <= 0
  ) {
    stop(
      mirna_id,
      " expression has an invalid range.",
      call. = FALSE
    )
  }

  bracket_rows <- lapply(
    seq_along(highlighted_contrasts),
    function(index) {
      contrast_name <- highlighted_contrasts[[index]]
      statistics_row <- selected_statistics[
        as.character(selected_statistics$contrast) ==
          contrast_name,
        ,
        drop = FALSE
      ]
      endpoints <- if (
        contrast_name == "MCAO3hr_vs_Sham"
      ) {
        c(1, 3)
      } else if (
        contrast_name == "MCAO3hr_vs_MCAO1hr"
      ) {
        c(2, 3)
      } else {
        c(1, 2)
      }
      bracket_height <- expression_range_local[[2L]] +
        (0.09 + (index - 1L) * 0.11) *
          expression_span_local
      data.frame(
        x_start = endpoints[[1L]],
        x_end = endpoints[[2L]],
        bracket_y = bracket_height,
        tip_y = bracket_height -
          0.025 * expression_span_local,
        label_y = bracket_height +
          0.038 * expression_span_local,
        label = paste0(
          contrast_labels[[contrast_name]],
          ": FDR = ",
          sprintf("%.3f", statistics_row$adj.P.Val)
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  bracket_data <- do.call(rbind, bracket_rows)
  plot_upper_local <- max(bracket_data$label_y) +
    0.07 * expression_span_local
  plot_lower_local <- expression_range_local[[1L]] -
    0.06 * expression_span_local
  subtitle_statistics <- selected_statistics[1L, , drop = FALSE]

  panel_plot <- ggplot2::ggplot(
    sample_expression_local,
    ggplot2::aes(
      x = treatment,
      y = expression,
      fill = treatment
    )
  ) +
    ggplot2::geom_violin(
      width = 0.72,
      scale = "width",
      trim = FALSE,
      alpha = 0.11,
      color = NA
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(
        width = 0.12,
        height = 0,
        seed = jitter_seed
      ),
      shape = 21,
      size = 3.0,
      stroke = 0.38,
      color = "white",
      alpha = 0.82
    ) +
    ggplot2::geom_line(
      data = adjusted_estimates_local,
      ggplot2::aes(
        x = treatment,
        y = adjusted_estimate,
        group = 1
      ),
      inherit.aes = FALSE,
      linewidth = 0.85,
      color = "#34424C"
    ) +
    ggplot2::geom_errorbar(
      data = adjusted_estimates_local,
      ggplot2::aes(
        x = treatment,
        ymin = confidence_low,
        ymax = confidence_high
      ),
      inherit.aes = FALSE,
      width = 0.13,
      linewidth = 0.85,
      color = "#202B33"
    ) +
    ggplot2::geom_point(
      data = adjusted_estimates_local,
      ggplot2::aes(
        x = treatment,
        y = adjusted_estimate,
        fill = treatment
      ),
      inherit.aes = FALSE,
      shape = 23,
      size = 5.1,
      stroke = 0.95,
      color = "#202B33"
    ) +
    ggplot2::scale_fill_manual(
      values = treatment_colors,
      guide = "none"
    ) +
    ggplot2::scale_x_discrete(
      labels = treatment_axis_labels_local
    ) +
    ggplot2::coord_cartesian(
      ylim = c(plot_lower_local, plot_upper_local),
      clip = "off"
    ) +
    ggplot2::labs(
      title = paste0(
        "(",
        panel_label,
        ") ",
        mirna_id,
        " expression"
      ),
      subtitle = paste0(
        contrast_labels[[
          as.character(subtitle_statistics$contrast)
        ]],
        "\nlog2FC = ",
        sprintf("%+.2f", subtitle_statistics$logFC),
        "; p = ",
        format_p_value(subtitle_statistics$P.Value),
        "; FDR = ",
        sprintf("%.3f", subtitle_statistics$adj.P.Val)
      ),
      x = NULL,
      y = expression(log[2] * "-normalized expression"),
      caption = paste0(
        "Circles: individual animals. Diamonds/error bars: adjusted ",
        "means and 95% CIs."
      )
    ) +
    ggplot2::theme_minimal(
      base_size = 15,
      base_family = "sans"
    ) +
    expression_plot$theme

  for (index in seq_len(nrow(bracket_data))) {
    bracket <- bracket_data[index, , drop = FALSE]
    panel_plot <- panel_plot +
      ggplot2::annotate(
        "segment",
        x = bracket$x_start,
        xend = bracket$x_end,
        y = bracket$bracket_y,
        yend = bracket$bracket_y,
        linewidth = 0.75,
        color = "#34424C"
      ) +
      ggplot2::annotate(
        "segment",
        x = c(bracket$x_start, bracket$x_end),
        xend = c(bracket$x_start, bracket$x_end),
        y = bracket$bracket_y,
        yend = bracket$tip_y,
        linewidth = 0.75,
        color = "#34424C"
      ) +
      ggplot2::annotate(
        "text",
        x = mean(c(bracket$x_start, bracket$x_end)),
        y = bracket$label_y,
        label = bracket$label,
        fontface = "bold",
        size = 4.8,
        color = "#283640"
      )
  }

  grDevices::png(
    filename = output_figure_file,
    width = 6.4,
    height = 6.7,
    units = "in",
    res = 400,
    bg = "white"
  )
  print(panel_plot)
  grDevices::dev.off()

  invisible(
    list(
      plot = panel_plot,
      statistics = selected_statistics,
      sample_expression = sample_expression_local,
      adjusted_estimates = adjusted_estimates_local
    )
  )
}

mirna_8101_panel <- make_additional_expression_panel(
  mirna_id = "miR-8101",
  panel_label = "C",
  highlighted_contrasts = c(
    "MCAO3hr_vs_Sham",
    "MCAO3hr_vs_MCAO1hr"
  ),
  output_figure_file = mirna_8101_figure_file,
  output_sample_file = mirna_8101_sample_file,
  output_adjusted_file = mirna_8101_adjusted_file,
  jitter_seed = 20260730
)

mirna_5130_panel <- make_additional_expression_panel(
  mirna_id = "miR-5130",
  panel_label = "D",
  highlighted_contrasts = c(
    "MCAO3hr_vs_Sham",
    "MCAO3hr_vs_MCAO1hr"
  ),
  output_figure_file = mirna_5130_figure_file,
  output_sample_file = mirna_5130_sample_file,
  output_adjusted_file = mirna_5130_adjusted_file,
  jitter_seed = 20260731
)

draw_main_figure_3 <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2L,
    ncol = 3L,
    heights = grid::unit(c(8.8, 6.7), "null"),
    widths = grid::unit(c(1, 1, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(
    coefficient_plot,
    vp = grid::viewport(
      layout.pos.row = 1L,
      layout.pos.col = 1:3
    )
  )
  print(
    expression_plot + ggplot2::labs(caption = NULL),
    vp = grid::viewport(
      layout.pos.row = 2L,
      layout.pos.col = 1L
    )
  )
  print(
    mirna_8101_panel$plot + ggplot2::labs(caption = NULL),
    vp = grid::viewport(
      layout.pos.row = 2L,
      layout.pos.col = 2L
    )
  )
  print(
    mirna_5130_panel$plot + ggplot2::labs(caption = NULL),
    vp = grid::viewport(
      layout.pos.row = 2L,
      layout.pos.col = 3L
    )
  )
  grid::popViewport()
}

grDevices::png(
  filename = main_figure_file,
  width = 17.5,
  height = 15.5,
  units = "in",
  res = 400,
  bg = "white"
)
draw_main_figure_3()
grDevices::dev.off()

message("Wrote Figure 3A: ", figure_file)
message("Wrote Figure 3A source data: ", source_data_file)
message("Wrote Figure 3B: ", expression_figure_file)
message("Wrote Figure 3B sample expression: ", sample_expression_file)
message("Wrote Figure 3B adjusted estimates: ", adjusted_estimate_file)
message("Wrote Figure 3C: ", mirna_8101_figure_file)
message("Wrote Figure 3C sample expression: ", mirna_8101_sample_file)
message("Wrote Figure 3C adjusted estimates: ", mirna_8101_adjusted_file)
message("Wrote Figure 3D: ", mirna_5130_figure_file)
message("Wrote Figure 3D sample expression: ", mirna_5130_sample_file)
message("Wrote Figure 3D adjusted estimates: ", mirna_5130_adjusted_file)
message("Wrote Main Figure 3: ", main_figure_file)
message(
  "Selected miRNAs: ",
  length(selected_feature_ids),
  "; residual df used for intervals: ",
  residual_degrees_freedom
)
