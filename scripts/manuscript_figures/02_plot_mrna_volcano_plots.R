#!/usr/bin/env Rscript

# Create mRNA volcano and adjusted expression plots for manuscript Figure 2.
#
# This script visualizes saved covariate-adjusted limma results. It does not
# rerun genome-wide differential expression. Model-adjusted group estimates
# for six selected genes are calculated from the saved expression matrix and
# the exact design used in the differential-expression analysis.
#
# Inputs:
#   results/mrna/differential_expression/treatment_sex_age_limma/
#     analysis_samples_mrna.csv
#     design_matrix_mrna.csv
#     transcript_cluster_MCAO1hr_vs_Sham_mrna.csv
#     transcript_cluster_MCAO3hr_vs_Sham_mrna.csv
#     transcript_cluster_MCAO3hr_vs_MCAO1hr_mrna.csv
#   results/mrna/analysis_ready/
#     expression_matrix_unique_gene_mapped_mrna.csv
#
# Outputs:
#   results/manuscript/figures/
#     main/figure2_mrna_differential_expression.png
#     components/figure2_mrna_volcano_plots.png
#     components/figure2_mrna_adjusted_expression_plots.png
#     source_data/figure2_mrna_volcano_summary.csv
#     source_data/figure2_mrna_volcano_high_confidence_genes.csv
#     source_data/figure2_mrna_adjusted_expression_estimates.csv
#
# Thresholds:
#   raw p-value < 0.05
#   BH FDR < 0.01
#   absolute log2 fold change >= 1
#
# Usage:
#   Rscript scripts/manuscript_figures/02_plot_mrna_volcano_plots.R

required_packages <- c("ggplot2", "ggrepel", "limma")
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
  "mrna",
  "differential_expression",
  "treatment_sex_age_limma"
)
figure_root <- file.path("results", "manuscript", "figures")
main_output_dir <- file.path(figure_root, "main")
component_output_dir <- file.path(figure_root, "components")
source_data_output_dir <- file.path(figure_root, "source_data")
output_directories <- c(
  main_output_dir,
  component_output_dir,
  source_data_output_dir
)
for (directory in output_directories) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

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
panel_labels <- c(
  MCAO1hr_vs_Sham = "A",
  MCAO3hr_vs_Sham = "B",
  MCAO3hr_vs_MCAO1hr = "C"
)

raw_p_cutoff <- 0.05
fdr_cutoff <- 0.01
absolute_logfc_cutoff <- 1
labels_per_contrast <- 18L

category_order <- c(
  "Other",
  "Downregulated",
  "Upregulated"
)
category_colors <- c(
  "Other" = "#C7CDD3",
  "Downregulated" = "#2166AC",
  "Upregulated" = "#D73027"
)

input_files <- stats::setNames(
  file.path(
    input_dir,
    paste0("transcript_cluster_", contrast_order, "_mrna.csv")
  ),
  contrast_order
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required mRNA differential-expression input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

plot_data <- list()
summary_rows <- list()
high_confidence_rows <- list()
label_rows <- list()

for (contrast_name in contrast_order) {
  result <- read.csv(
    input_files[[contrast_name]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c(
    "transcript_cluster_id",
    "SYMBOL",
    "logFC",
    "P.Value",
    "adj.P.Val"
  )
  missing_columns <- setdiff(required_columns, names(result))
  if (length(missing_columns) > 0L) {
    stop(
      basename(input_files[[contrast_name]]),
      " lacks required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(result$transcript_cluster_id)) {
    stop(
      basename(input_files[[contrast_name]]),
      " contains duplicated transcript_cluster_id values.",
      call. = FALSE
    )
  }
  if (
    any(!is.finite(result$logFC)) ||
    any(!is.finite(result$P.Value)) ||
    any(!is.finite(result$adj.P.Val))
  ) {
    stop(
      basename(input_files[[contrast_name]]),
      " contains non-finite differential-expression statistics.",
      call. = FALSE
    )
  }
  if (nrow(result) != 20222L) {
    stop(
      basename(input_files[[contrast_name]]),
      " does not contain the expected 20,222 tested features.",
      call. = FALSE
    )
  }

  result$contrast <- contrast_name
  result$minus_log10_p <- -log10(
    pmax(result$P.Value, .Machine$double.xmin)
  )
  result$category <- "Other"
  result$category[
    result$P.Value < raw_p_cutoff &
      result$logFC <= -absolute_logfc_cutoff
  ] <- "Downregulated"
  result$category[
    result$P.Value < raw_p_cutoff &
      result$logFC >= absolute_logfc_cutoff
  ] <- "Upregulated"
  result$category <- factor(
    result$category,
    levels = category_order
  )
  result <- result[
    order(result$category),
    ,
    drop = FALSE
  ]
  plot_data[[contrast_name]] <- result

  combined_threshold <- (
    result$adj.P.Val < fdr_cutoff &
      abs(result$logFC) >= absolute_logfc_cutoff
  )
  high_confidence <- result[combined_threshold, , drop = FALSE]
  high_confidence$gene_label <- ifelse(
    is.na(high_confidence$SYMBOL) |
      trimws(high_confidence$SYMBOL) == "",
    high_confidence$transcript_cluster_id,
    high_confidence$SYMBOL
  )
  high_confidence_rows[[contrast_name]] <- high_confidence

  label_candidates <- result[
    result$P.Value < raw_p_cutoff &
      abs(result$logFC) >= absolute_logfc_cutoff,
    ,
    drop = FALSE
  ]
  label_candidates <- label_candidates[
    order(label_candidates$P.Value, -abs(label_candidates$logFC)),
    ,
    drop = FALSE
  ]
  label_candidates <- utils::head(
    label_candidates,
    labels_per_contrast
  )
  label_candidates$gene_label <- ifelse(
    is.na(label_candidates$SYMBOL) |
      trimws(label_candidates$SYMBOL) == "",
    label_candidates$transcript_cluster_id,
    label_candidates$SYMBOL
  )
  label_rows[[contrast_name]] <- label_candidates

  summary_rows[[contrast_name]] <- data.frame(
    contrast = contrast_name,
    tested_features = nrow(result),
    raw_p_lt_0_05 = sum(result$P.Value < raw_p_cutoff),
    fdr_lt_0_01 = sum(result$adj.P.Val < fdr_cutoff),
    raw_p_and_abs_logfc_ge_1 = sum(
      result$P.Value < raw_p_cutoff &
        abs(result$logFC) >= absolute_logfc_cutoff
    ),
    fdr_and_abs_logfc_ge_1 = sum(combined_threshold),
    fdr_and_logfc_le_minus_1 = sum(
      result$adj.P.Val < fdr_cutoff &
        result$logFC <= -absolute_logfc_cutoff
    ),
    fdr_and_logfc_ge_1 = sum(
      result$adj.P.Val < fdr_cutoff &
        result$logFC >= absolute_logfc_cutoff
    ),
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL
expected_combined_counts <- c(0L, 0L, 0L)
if (
  !identical(
    as.integer(summary_table$fdr_and_abs_logfc_ge_1),
    expected_combined_counts
  )
) {
  stop(
    paste0(
      "The number of FDR-supported genes with at least two-fold change ",
      "differs from the expected saved results."
    ),
    call. = FALSE
  )
}

all_x <- unlist(lapply(plot_data, `[[`, "logFC"), use.names = FALSE)
all_y <- unlist(
  lapply(plot_data, `[[`, "minus_log10_p"),
  use.names = FALSE
)
x_limit <- max(2.4, ceiling(10 * max(abs(all_x))) / 10)
y_limit <- ceiling(10 * max(all_y + 0.35)) / 10

volcano_theme <- ggplot2::theme_minimal(
  base_size = 15,
  base_family = "sans"
) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 17.5,
      color = "#18212A"
    ),
    plot.subtitle = ggplot2::element_text(
      face = "bold",
      size = 13.2,
      color = "#52616D",
      margin = ggplot2::margin(b = 7)
    ),
    axis.title = ggplot2::element_text(
      face = "bold",
      size = 15,
      color = "#28333D"
    ),
    axis.text = ggplot2::element_text(
      face = "bold",
      size = 13.2,
      color = "#34414D"
    ),
    panel.grid.major = ggplot2::element_line(
      color = "#E6E6E6",
      linewidth = 0.32
    ),
    panel.grid.minor = ggplot2::element_line(
      color = "#F2F2F2",
      linewidth = 0.18
    ),
    panel.border = ggplot2::element_rect(
      color = "#56636D",
      fill = NA,
      linewidth = 0.55
    ),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(
      face = "bold",
      size = 13.2,
      color = "#18212A"
    ),
    plot.margin = ggplot2::margin(8, 9, 7, 9)
  )

make_volcano_plot <- function(contrast_name) {
  result <- plot_data[[contrast_name]]
  high_confidence <- high_confidence_rows[[contrast_name]]
  label_candidates <- label_rows[[contrast_name]]
  combined_count <- nrow(high_confidence)

  plot <- ggplot2::ggplot(
    result,
    ggplot2::aes(
      x = logFC,
      y = minus_log10_p,
      color = category
    )
  ) +
    ggplot2::geom_hline(
      yintercept = -log10(raw_p_cutoff),
      linetype = 2,
      linewidth = 0.6,
      color = "#58656F"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = 3,
      linewidth = 0.72,
      color = "#6A3D9A"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-absolute_logfc_cutoff, absolute_logfc_cutoff),
      linetype = 2,
      linewidth = 0.6,
      color = "#58656F"
    ) +
    ggplot2::geom_point(
      size = 1.75,
      alpha = 0.66,
      stroke = 0
    ) +
    ggplot2::geom_point(
      data = high_confidence,
      size = 3.15,
      alpha = 0.98,
      stroke = 0
    ) +
    ggplot2::scale_color_manual(
      values = category_colors,
      breaks = category_order,
      drop = FALSE
    ) +
    ggplot2::coord_cartesian(
      xlim = c(-x_limit, x_limit),
      ylim = c(0, y_limit),
      clip = "off"
    ) +
    ggplot2::labs(
      title = paste0(
        "(",
        panel_labels[[contrast_name]],
        ")  ",
        contrast_labels[[contrast_name]]
      ),
      subtitle = paste0(
        "FDR < 0.01 and |log2FC| ≥ 1: n = ",
        combined_count
      ),
      x = expression(log[2] * " fold change"),
      y = expression(-log[10] * "(p-value)")
    ) +
    volcano_theme

  if (nrow(label_candidates) > 0L) {
    plot <- plot +
      ggrepel::geom_text_repel(
        data = label_candidates,
        ggplot2::aes(label = gene_label),
        size = 3.4,
        fontface = "bold",
        color = "#18212A",
        seed = 20260729,
        box.padding = 0.65,
        point.padding = 0.30,
        force = 2.4,
        force_pull = 0.45,
        max.time = 5,
        max.iter = 50000,
        min.segment.length = 0,
        segment.color = "#7B8790",
        segment.size = 0.42,
        max.overlaps = Inf,
        show.legend = FALSE
      )
  }
  plot
}

volcano_plots <- lapply(contrast_order, make_volcano_plot)

draw_volcano_composite <- function() {
  legend_source <- ggplot2::ggplotGrob(
    volcano_plots[[3]] + ggplot2::theme(legend.position = "bottom")
  )
  legend_index <- which(
    legend_source$layout$name == "guide-box-bottom"
  )
  if (length(legend_index) != 1L) {
    stop("Could not extract the shared volcano-plot legend.", call. = FALSE)
  }
  shared_legend <- legend_source$grobs[[legend_index]]
  plots_without_legends <- lapply(
    volcano_plots,
    function(plot) plot + ggplot2::theme(legend.position = "none")
  )

  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2L,
    ncol = 3L,
    widths = grid::unit(c(1, 1, 1), "null"),
    heights = grid::unit(c(1, 0.13), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  for (column in seq_along(plots_without_legends)) {
    print(
      plots_without_legends[[column]],
      vp = grid::viewport(
        layout.pos.row = 1L,
        layout.pos.col = column
      )
    )
  }
  grid::pushViewport(
    grid::viewport(layout.pos.row = 2L, layout.pos.col = 1:3)
  )
  grid::grid.draw(shared_legend)
  grid::popViewport()
  grid::popViewport()
}

figure_file <- file.path(
  component_output_dir,
  "figure2_mrna_volcano_plots.png"
)
grDevices::png(
  filename = figure_file,
  width = 16.2,
  height = 7.2,
  units = "in",
  res = 400,
  bg = "white"
)
draw_volcano_composite()
grDevices::dev.off()

write.csv(
  summary_table,
  file.path(
    source_data_output_dir,
    "figure2_mrna_volcano_summary.csv"
  ),
  row.names = FALSE
)

high_confidence_table <- do.call(
  rbind,
  lapply(
    high_confidence_rows,
    function(result) {
      result[
        ,
        c(
          "contrast",
          "transcript_cluster_id",
          "SYMBOL",
          "logFC",
          "P.Value",
          "adj.P.Val"
        )
      ]
    }
  )
)
rownames(high_confidence_table) <- NULL
write.csv(
  high_confidence_table,
  file.path(
    source_data_output_dir,
    "figure2_mrna_volcano_high_confidence_genes.csv"
  ),
  row.names = FALSE
)

selected_gene_specs <- data.frame(
  panel = LETTERS[4:9],
  gene = c(
    "Ptprcap",
    "Srsf10",
    "Dusp5",
    "Il1r2",
    "Pram1",
    "Mmp16"
  ),
  highlighted_contrast = c(
    "MCAO3hr_vs_MCAO1hr",
    "MCAO3hr_vs_MCAO1hr",
    "MCAO3hr_vs_MCAO1hr",
    "MCAO3hr_vs_Sham",
    "MCAO3hr_vs_Sham",
    "MCAO3hr_vs_MCAO1hr"
  ),
  stringsAsFactors = FALSE
)

selected_result_rows <- lapply(
  seq_len(nrow(selected_gene_specs)),
  function(index) {
    specification <- selected_gene_specs[index, , drop = FALSE]
    result <- plot_data[[specification$highlighted_contrast]]
    selected <- result[
      result$SYMBOL == specification$gene,
      ,
      drop = FALSE
    ]
    if (nrow(selected) != 1L) {
      stop(
        "Expected exactly one measured transcript cluster for ",
        specification$gene,
        " in ",
        specification$highlighted_contrast,
        ".",
        call. = FALSE
      )
    }
    selected
  }
)
selected_statistics <- do.call(rbind, selected_result_rows)
rownames(selected_statistics) <- NULL
selected_gene_specs$transcript_cluster_id <-
  selected_statistics$transcript_cluster_id
selected_gene_specs$logFC <- selected_statistics$logFC
selected_gene_specs$P.Value <- selected_statistics$P.Value
selected_gene_specs$adj.P.Val <- selected_statistics$adj.P.Val

expression_file <- file.path(
  "results",
  "mrna",
  "analysis_ready",
  "expression_matrix_unique_gene_mapped_mrna.csv"
)
analysis_sample_file <- file.path(
  input_dir,
  "analysis_samples_mrna.csv"
)
design_file <- file.path(input_dir, "design_matrix_mrna.csv")
expression_inputs <- c(
  expression_file,
  analysis_sample_file,
  design_file
)
missing_expression_inputs <- expression_inputs[
  !file.exists(expression_inputs)
]
if (length(missing_expression_inputs) > 0L) {
  stop(
    "Missing required adjusted-expression input(s): ",
    paste(missing_expression_inputs, collapse = ", "),
    call. = FALSE
  )
}

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
  "age_group"
)
missing_sample_columns <- setdiff(
  required_sample_columns,
  names(analysis_samples)
)
if (length(missing_sample_columns) > 0L) {
  stop(
    "Saved analysis samples lack required column(s): ",
    paste(missing_sample_columns, collapse = ", "),
    call. = FALSE
  )
}
if (
  nrow(analysis_samples) != 43L ||
    anyDuplicated(analysis_samples$file_name) ||
    anyDuplicated(analysis_samples$animal_id)
) {
  stop(
    "Expected 43 unique paired animals in the saved mRNA analysis samples.",
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
analysis_samples$age_group <- factor(
  analysis_samples$age_group,
  levels = c("Old", "Young")
)
if (
  anyNA(analysis_samples$treatment) ||
    anyNA(analysis_samples$sex) ||
    anyNA(analysis_samples$age_group)
) {
  stop(
    "Treatment, sex, or age is missing from the saved mRNA analysis samples.",
    call. = FALSE
  )
}

saved_design <- as.matrix(
  read.csv(
    design_file,
    row.names = 1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
)
storage.mode(saved_design) <- "numeric"
rebuilt_design <- stats::model.matrix(
  ~ 0 + treatment + sex + age_group,
  data = analysis_samples
)
colnames(rebuilt_design) <- sub(
  "^treatment",
  "",
  colnames(rebuilt_design)
)
if (
  !identical(dim(saved_design), dim(rebuilt_design)) ||
    !identical(colnames(saved_design), colnames(rebuilt_design)) ||
    !isTRUE(
      all.equal(
        as.numeric(saved_design),
        as.numeric(rebuilt_design),
        tolerance = 1e-12
      )
    )
) {
  stop(
    "The saved design matrix does not match treatment + sex + age.",
    call. = FALSE
  )
}

expression_table <- read.csv(
  expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (
  !"transcript_cluster_id" %in% names(expression_table) ||
    !all(analysis_samples$file_name %in% names(expression_table))
) {
  stop(
    "The analysis-ready mRNA matrix lacks required identifiers or samples.",
    call. = FALSE
  )
}
selected_expression_index <- match(
  selected_gene_specs$transcript_cluster_id,
  expression_table$transcript_cluster_id
)
if (anyNA(selected_expression_index)) {
  stop(
    "At least one selected transcript cluster is absent from the expression matrix.",
    call. = FALSE
  )
}
selected_expression <- as.matrix(
  expression_table[
    selected_expression_index,
    analysis_samples$file_name,
    drop = FALSE
  ]
)
storage.mode(selected_expression) <- "numeric"
rownames(selected_expression) <- selected_gene_specs$gene
if (any(!is.finite(selected_expression))) {
  stop(
    "Selected-gene expression contains non-finite values.",
    call. = FALSE
  )
}

selected_fit <- limma::lmFit(selected_expression, saved_design)
covariate_columns <- c("sexMale", "age_groupYoung")
if (!all(covariate_columns %in% colnames(saved_design))) {
  stop(
    "The saved model lacks the expected sex and age covariates.",
    call. = FALSE
  )
}
covariate_means <- colMeans(
  saved_design[, covariate_columns, drop = FALSE]
)

adjusted_rows <- list()
raw_expression_rows <- list()
for (index in seq_len(nrow(selected_gene_specs))) {
  specification <- selected_gene_specs[index, , drop = FALSE]
  gene <- specification$gene
  coefficient <- selected_fit$coefficients[gene, ]

  if (specification$highlighted_contrast == "MCAO3hr_vs_Sham") {
    fitted_logfc <- coefficient[["MCAO3hr"]] - coefficient[["Sham"]]
  } else if (
    specification$highlighted_contrast == "MCAO3hr_vs_MCAO1hr"
  ) {
    fitted_logfc <- coefficient[["MCAO3hr"]] -
      coefficient[["MCAO1hr"]]
  } else {
    stop("Unsupported highlighted contrast.", call. = FALSE)
  }
  if (
    !isTRUE(
      all.equal(
        unname(fitted_logfc),
        specification$logFC,
        tolerance = 1e-8
      )
    )
  ) {
    stop(
      "Selected-gene model coefficient does not reproduce the saved logFC for ",
      gene,
      ".",
      call. = FALSE
    )
  }

  raw_expression_rows[[gene]] <- data.frame(
    gene = gene,
    panel = specification$panel,
    animal_id = analysis_samples$animal_id,
    treatment = analysis_samples$treatment,
    expression = as.numeric(selected_expression[gene, ]),
    stringsAsFactors = FALSE
  )

  for (treatment_name in treatment_order) {
    estimate_vector <- stats::setNames(
      rep(0, ncol(saved_design)),
      colnames(saved_design)
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
    standard_error <- selected_fit$sigma[[gene]] *
      unscaled_standard_error
    degrees_freedom <- if (length(selected_fit$df.residual) == 1L) {
      selected_fit$df.residual[[1L]]
    } else {
      selected_fit$df.residual[[index]]
    }
    critical_value <- stats::qt(0.975, df = degrees_freedom)

    adjusted_rows[[paste(gene, treatment_name, sep = "_")]] <- data.frame(
      gene = gene,
      panel = specification$panel,
      treatment = treatment_name,
      adjusted_estimate = estimate,
      standard_error = standard_error,
      confidence_low = estimate - critical_value * standard_error,
      confidence_high = estimate + critical_value * standard_error,
      highlighted_contrast = specification$highlighted_contrast,
      logFC = specification$logFC,
      P.Value = specification$P.Value,
      adj.P.Val = specification$adj.P.Val,
      stringsAsFactors = FALSE
    )
  }
}

raw_expression_data <- do.call(rbind, raw_expression_rows)
rownames(raw_expression_data) <- NULL
raw_expression_data$treatment <- factor(
  raw_expression_data$treatment,
  levels = treatment_order
)
adjusted_estimates <- do.call(rbind, adjusted_rows)
rownames(adjusted_estimates) <- NULL
adjusted_estimates$treatment <- factor(
  adjusted_estimates$treatment,
  levels = treatment_order
)

adjusted_estimate_file <- file.path(
  source_data_output_dir,
  "figure2_mrna_adjusted_expression_estimates.csv"
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
highlighted_contrast_labels <- c(
  MCAO3hr_vs_Sham = "MCAO3hr vs Sham",
  MCAO3hr_vs_MCAO1hr = "MCAO3hr vs MCAO1hr"
)

format_fdr <- function(value) {
  if (value < 0.001) {
    format(value, scientific = TRUE, digits = 2)
  } else {
    sprintf("%.3f", value)
  }
}

expression_theme <- ggplot2::theme_minimal(
  base_size = 15,
  base_family = "sans"
) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 17,
      color = "#18212A"
    ),
    plot.subtitle = ggplot2::element_text(
      face = "bold",
      size = 13,
      color = "#52616D",
      margin = ggplot2::margin(b = 6)
    ),
    axis.title = ggplot2::element_text(
      face = "bold",
      size = 15,
      color = "#28333D"
    ),
    axis.text = ggplot2::element_text(
      face = "bold",
      size = 13.2,
      color = "#34414D"
    ),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E7EBEE",
      linewidth = 0.34
    ),
    panel.border = ggplot2::element_rect(
      color = "#56636D",
      fill = NA,
      linewidth = 0.55
    ),
    legend.position = "none",
    plot.margin = ggplot2::margin(8, 10, 8, 10)
  )

make_expression_plot <- function(index) {
  specification <- selected_gene_specs[index, , drop = FALSE]
  gene <- specification$gene
  raw_gene <- raw_expression_data[
    raw_expression_data$gene == gene,
    ,
    drop = FALSE
  ]
  adjusted_gene <- adjusted_estimates[
    adjusted_estimates$gene == gene,
    ,
    drop = FALSE
  ]
  subtitle <- paste0(
    highlighted_contrast_labels[[specification$highlighted_contrast]],
    ": log2FC = ",
    sprintf("%+.2f", specification$logFC),
    "; FDR = ",
    format_fdr(specification$adj.P.Val)
  )

  ggplot2::ggplot(
    raw_gene,
    ggplot2::aes(
      x = treatment,
      y = expression,
      color = treatment
    )
  ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(
        width = 0.13,
        height = 0,
        seed = 20260729 + index
      ),
      size = 2.35,
      alpha = 0.68,
      stroke = 0
    ) +
    ggplot2::geom_line(
      data = adjusted_gene,
      ggplot2::aes(
        x = treatment,
        y = adjusted_estimate,
        group = 1
      ),
      inherit.aes = FALSE,
      linewidth = 0.72,
      color = "#3C4852"
    ) +
    ggplot2::geom_errorbar(
      data = adjusted_gene,
      ggplot2::aes(
        x = treatment,
        ymin = confidence_low,
        ymax = confidence_high
      ),
      inherit.aes = FALSE,
      width = 0.12,
      linewidth = 0.78,
      color = "#202A33"
    ) +
    ggplot2::geom_point(
      data = adjusted_gene,
      ggplot2::aes(
        x = treatment,
        y = adjusted_estimate,
        fill = treatment
      ),
      inherit.aes = FALSE,
      shape = 23,
      size = 4.1,
      stroke = 0.85,
      color = "#202A33"
    ) +
    ggplot2::scale_color_manual(values = treatment_colors) +
    ggplot2::scale_fill_manual(values = treatment_colors) +
    ggplot2::labs(
      title = paste0(
        "(",
        specification$panel,
        ")  ",
        gene
      ),
      subtitle = subtitle,
      x = NULL,
      y = expression(log[2] * "-normalized expression")
    ) +
    expression_theme
}

expression_plots <- lapply(
  seq_len(nrow(selected_gene_specs)),
  make_expression_plot
)

draw_expression_composite <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2L,
    ncol = 3L,
    widths = grid::unit(c(1, 1, 1), "null"),
    heights = grid::unit(c(1, 1), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  for (index in seq_along(expression_plots)) {
    row <- ceiling(index / 3)
    column <- ((index - 1L) %% 3L) + 1L
    print(
      expression_plots[[index]],
      vp = grid::viewport(
        layout.pos.row = row,
        layout.pos.col = column
      )
    )
  }
  grid::popViewport()
}

expression_figure_file <- file.path(
  component_output_dir,
  "figure2_mrna_adjusted_expression_plots.png"
)
grDevices::png(
  filename = expression_figure_file,
  width = 16.2,
  height = 9.4,
  units = "in",
  res = 400,
  bg = "white"
)
draw_expression_composite()
grDevices::dev.off()

extract_volcano_legend <- function() {
  legend_source <- ggplot2::ggplotGrob(
    volcano_plots[[3]] + ggplot2::theme(legend.position = "bottom")
  )
  legend_index <- which(
    legend_source$layout$name == "guide-box-bottom"
  )
  if (length(legend_index) != 1L) {
    stop("Could not extract the volcano-plot legend.", call. = FALSE)
  }
  legend_source$grobs[[legend_index]]
}

draw_full_figure2 <- function() {
  shared_legend <- extract_volcano_legend()
  volcano_without_legends <- lapply(
    volcano_plots,
    function(plot) plot + ggplot2::theme(legend.position = "none")
  )
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 4L,
    ncol = 3L,
    widths = grid::unit(c(1, 1, 1), "null"),
    heights = grid::unit(c(1.18, 0.12, 0.93, 0.93), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  for (column in seq_along(volcano_without_legends)) {
    print(
      volcano_without_legends[[column]],
      vp = grid::viewport(
        layout.pos.row = 1L,
        layout.pos.col = column
      )
    )
  }
  grid::pushViewport(
    grid::viewport(layout.pos.row = 2L, layout.pos.col = 1:3)
  )
  grid::grid.draw(shared_legend)
  grid::popViewport()
  for (index in seq_along(expression_plots)) {
    row <- ceiling(index / 3) + 2L
    column <- ((index - 1L) %% 3L) + 1L
    print(
      expression_plots[[index]],
      vp = grid::viewport(
        layout.pos.row = row,
        layout.pos.col = column
      )
    )
  }
  grid::popViewport()
}

combined_figure_file <- file.path(
  main_output_dir,
  "figure2_mrna_differential_expression.png"
)
grDevices::png(
  filename = combined_figure_file,
  width = 16.2,
  height = 15.8,
  units = "in",
  res = 400,
  bg = "white"
)
draw_full_figure2()
grDevices::dev.off()

message("Wrote mRNA volcano plots: ", figure_file)
message(
  "Wrote adjusted expression plots: ",
  expression_figure_file
)
message("Wrote complete Figure 2: ", combined_figure_file)
