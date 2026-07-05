#!/usr/bin/env Rscript

# Create the paired-cohort PCA panels for manuscript Figure 1.
#
# This script performs visualization only. It does not rerun normalization,
# differential expression, association testing, or pathway analysis.
#
# Inputs:
#   results/multiomics/sample_manifest/
#     paired_manifest_mirna_mrna.csv
#   results/mrna/analysis_ready/
#     expression_matrix_unique_gene_mapped_mrna.csv
#   results/mirna/expression/rma_normalized_mirna/annotation/
#     mouse_mature_mirna_expression_mirna.csv
#   results/mirna/differential_expression/rma_normalized_mirna/
#     treatment_sex_age_limma/dabg_detection_summary_mirna.csv
#
# Outputs:
#   results/manuscript/figures/
#     main/figure1_experimental_design_and_paired_pca.png
#     supplementary/figure_s1_pca_stratified_by_sex.png
#     supplementary/figure_s2_pca_stratified_by_age.png
#     components/figure1_paired_data_overview.png
#     source_data/figure1_pca_scores.csv
#     source_data/figure1_summary.txt
#
# Usage:
#   Rscript scripts/manuscript_figures/01_plot_paired_data_overview.R

required_packages <- "ggplot2"
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

manifest_file <- file.path(
  "results", "multiomics", "sample_manifest",
  "paired_manifest_mirna_mrna.csv"
)
mrna_expression_file <- file.path(
  "results", "mrna", "analysis_ready",
  "expression_matrix_unique_gene_mapped_mrna.csv"
)
mirna_expression_file <- file.path(
  "results", "mirna", "expression", "rma_normalized_mirna",
  "annotation", "mouse_mature_mirna_expression_mirna.csv"
)
mirna_filter_file <- file.path(
  "results", "mirna", "differential_expression",
  "rma_normalized_mirna", "treatment_sex_age_limma",
  "dabg_detection_summary_mirna.csv"
)
figure_root <- file.path("results", "manuscript", "figures")
main_output_dir <- file.path(figure_root, "main")
supplementary_output_dir <- file.path(
  figure_root,
  "supplementary"
)
component_output_dir <- file.path(figure_root, "components")
source_data_output_dir <- file.path(figure_root, "source_data")

required_inputs <- c(
  manifest_file,
  mrna_expression_file,
  mirna_expression_file,
  mirna_filter_file
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}
output_directories <- c(
  main_output_dir,
  supplementary_output_dir,
  component_output_dir,
  source_data_output_dir
)
for (directory in output_directories) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

treatment_order <- c("Sham", "MCAO1hr", "MCAO3hr")
treatment_colors <- c(
  Sham = "#13D82F",
  MCAO1hr = "#1236D9",
  MCAO3hr = "#E5162E"
)
treatment_fill_colors <- c(
  Sham = "#78E88D",
  MCAO1hr = "#65B5FF",
  MCAO3hr = "#FF8F9D"
)
sex_order <- c("Female", "Male")
age_order <- c("Young", "Old")

publication_theme <- function(base_size = 10.5) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = base_size + 1.5, color = "#18212A"
      ),
      plot.subtitle = ggplot2::element_text(
        size = base_size - 0.5, color = "#4C5965", margin = ggplot2::margin(b = 8)
      ),
      axis.title = ggplot2::element_text(face = "bold", color = "#28333D"),
      axis.text = ggplot2::element_text(
        face = "bold",
        color = "#34414D"
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = base_size + 1.5,
        face = "bold",
        color = "#18212A"
      ),
      legend.key.width = grid::unit(1.2, "lines"),
      panel.background = ggplot2::element_rect(
        fill = "#FFFFFF",
        color = NA
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
        color = "#48545F", fill = NA, linewidth = 0.55
      ),
      plot.margin = ggplot2::margin(10, 12, 8, 10)
    )
}

draw_side_by_side <- function(left_plot, right_plot) {
  legend_source <- ggplot2::ggplotGrob(
    right_plot + ggplot2::theme(legend.position = "bottom")
  )
  legend_index <- which(
    legend_source$layout$name == "guide-box-bottom"
  )
  if (length(legend_index) != 1L) {
    stop("Could not extract the shared PCA legend.", call. = FALSE)
  }
  shared_legend <- legend_source$grobs[[legend_index]]
  left_plot <- left_plot + ggplot2::theme(legend.position = "none")
  right_plot <- right_plot + ggplot2::theme(legend.position = "none")

  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2L,
    ncol = 2L,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit(c(1, 0.10), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(
    left_plot,
    vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L)
  )
  print(
    right_plot,
    vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 2L)
  )
  grid::pushViewport(
    grid::viewport(layout.pos.row = 2L, layout.pos.col = 1:2)
  )
  grid::grid.draw(shared_legend)
  grid::popViewport()
  grid::popViewport()
}

save_composite_png <- function(
  left_plot,
  right_plot,
  stem,
  width,
  height,
  output_directory
) {
  grDevices::png(
    filename = file.path(
      output_directory,
      paste0(stem, ".png")
    ),
    width = width,
    height = height,
    units = "in",
    res = 400,
    bg = "white"
  )
  draw_side_by_side(left_plot, right_plot)
  grDevices::dev.off()
}

draw_stacked <- function(top_plot, bottom_plot) {
  legend_source <- ggplot2::ggplotGrob(
    top_plot + ggplot2::theme(legend.position = "bottom")
  )
  legend_index <- which(
    legend_source$layout$name == "guide-box-bottom"
  )
  if (length(legend_index) != 1L) {
    stop("Could not extract the shared PCA legend.", call. = FALSE)
  }
  shared_legend <- legend_source$grobs[[legend_index]]
  top_plot <- top_plot + ggplot2::theme(legend.position = "none")
  bottom_plot <- bottom_plot + ggplot2::theme(legend.position = "none")

  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 3L,
    ncol = 1L,
    heights = grid::unit(c(1, 1, 0.10), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(
    top_plot,
    vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L)
  )
  print(
    bottom_plot,
    vp = grid::viewport(layout.pos.row = 2L, layout.pos.col = 1L)
  )
  grid::pushViewport(
    grid::viewport(layout.pos.row = 3L, layout.pos.col = 1L)
  )
  grid::grid.draw(shared_legend)
  grid::popViewport()
  grid::popViewport()
}

save_stacked_png <- function(
  top_plot,
  bottom_plot,
  stem,
  width,
  height,
  output_directory
) {
  grDevices::png(
    filename = file.path(
      output_directory,
      paste0(stem, ".png")
    ),
    width = width,
    height = height,
    units = "in",
    res = 400,
    bg = "white"
  )
  draw_stacked(top_plot, bottom_plot)
  grDevices::dev.off()
}

draw_four_panel <- function(
  top_left_plot,
  top_right_plot,
  bottom_left_plot,
  bottom_right_plot
) {
  legend_source <- ggplot2::ggplotGrob(
    bottom_right_plot + ggplot2::theme(legend.position = "bottom")
  )
  legend_index <- which(
    legend_source$layout$name == "guide-box-bottom"
  )
  if (length(legend_index) != 1L) {
    stop("Could not extract the shared Figure 1 legend.", call. = FALSE)
  }
  shared_legend <- legend_source$grobs[[legend_index]]
  plots <- lapply(
    list(
      top_left_plot,
      top_right_plot,
      bottom_left_plot,
      bottom_right_plot
    ),
    function(plot) plot + ggplot2::theme(legend.position = "none")
  )

  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 3L,
    ncol = 2L,
    widths = grid::unit(c(1, 1), "null"),
    heights = grid::unit(c(0.58, 1, 0.09), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(
    plots[[1]],
    vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L)
  )
  print(
    plots[[2]],
    vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 2L)
  )
  print(
    plots[[3]],
    vp = grid::viewport(layout.pos.row = 2L, layout.pos.col = 1L)
  )
  print(
    plots[[4]],
    vp = grid::viewport(layout.pos.row = 2L, layout.pos.col = 2L)
  )
  grid::pushViewport(
    grid::viewport(layout.pos.row = 3L, layout.pos.col = 1:2)
  )
  grid::grid.draw(shared_legend)
  grid::popViewport()
  grid::popViewport()
}

save_four_panel_png <- function(
  top_left_plot,
  top_right_plot,
  bottom_left_plot,
  bottom_right_plot,
  stem,
  width,
  height,
  output_directory
) {
  grDevices::png(
    filename = file.path(
      output_directory,
      paste0(stem, ".png")
    ),
    width = width,
    height = height,
    units = "in",
    res = 400,
    bg = "white"
  )
  draw_four_panel(
    top_left_plot,
    top_right_plot,
    bottom_left_plot,
    bottom_right_plot
  )
  grDevices::dev.off()
}

manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)
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
if (nrow(manifest) != 43L || anyDuplicated(manifest$animal_id)) {
  stop(
    "Expected 43 unique paired animals in the manifest.",
    call. = FALSE
  )
}
manifest$treatment <- factor(
  manifest$treatment,
  levels = treatment_order
)
manifest$sex <- factor(
  manifest$sex,
  levels = sex_order
)
manifest$age_group <- factor(
  manifest$age_group,
  levels = age_order
)
if (anyNA(manifest$sex) || anyNA(manifest$age_group)) {
  stop(
    "Paired manifest contains missing or unexpected sex or age-group labels.",
    call. = FALSE
  )
}
observed_counts <- table(manifest$treatment)
expected_counts <- c(Sham = 16L, MCAO1hr = 16L, MCAO3hr = 11L)
if (!identical(as.integer(observed_counts), as.integer(expected_counts))) {
  stop(
    "Unexpected paired-cohort treatment counts.",
    call. = FALSE
  )
}

mrna_expression <- read.csv(
  mrna_expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!"transcript_cluster_id" %in% names(mrna_expression)) {
  stop("mRNA expression lacks transcript_cluster_id.", call. = FALSE)
}
mrna_columns <- match(manifest$mrna_file_name, names(mrna_expression))
if (anyNA(mrna_columns)) {
  stop(
    "Could not match every paired animal to the mRNA expression matrix.",
    call. = FALSE
  )
}
mrna_matrix <- as.matrix(mrna_expression[, mrna_columns, drop = FALSE])
storage.mode(mrna_matrix) <- "double"
if (any(!is.finite(mrna_matrix))) {
  stop("mRNA expression contains non-finite values.", call. = FALSE)
}

mirna_expression <- read.csv(
  mirna_expression_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
mirna_filter <- read.csv(
  mirna_filter_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_mirna_filter_columns <- c(
  "ProbeSetName", "dabg_light_filter_keep"
)
if (!all(required_mirna_filter_columns %in% names(mirna_filter))) {
  stop("miRNA filter table has an unexpected schema.", call. = FALSE)
}
kept_mirna_ids <- as.character(
  mirna_filter$ProbeSetName[mirna_filter$dabg_light_filter_keep]
)
mirna_expression <- mirna_expression[
  as.character(mirna_expression$ProbeSetName) %in% kept_mirna_ids,
  ,
  drop = FALSE
]
if (nrow(mirna_expression) != 467L) {
  stop(
    "Expected 467 reliably detected mature mouse miRNA probesets.",
    call. = FALSE
  )
}
mirna_columns <- match(manifest$mirna_file_name, names(mirna_expression))
if (anyNA(mirna_columns)) {
  stop(
    "Could not match every paired animal to the miRNA expression matrix.",
    call. = FALSE
  )
}
mirna_matrix <- as.matrix(mirna_expression[, mirna_columns, drop = FALSE])
storage.mode(mirna_matrix) <- "double"
if (any(!is.finite(mirna_matrix))) {
  stop("miRNA expression contains non-finite values.", call. = FALSE)
}

compute_pca_scores <- function(expression_matrix, layer) {
  pca <- stats::prcomp(
    t(expression_matrix),
    center = TRUE,
    scale. = FALSE
  )
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  scores <- data.frame(
    layer = layer,
    animal_id = manifest$animal_id,
    treatment = manifest$treatment,
    sex = manifest$sex,
    age_group = manifest$age_group,
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC1_variance_percent = variance[[1]],
    PC2_variance_percent = variance[[2]],
    stringsAsFactors = FALSE
  )
  list(scores = scores, variance = variance)
}

add_point_depth_coordinates <- function(scores) {
  x_span <- diff(range(scores$PC1, finite = TRUE))
  y_span <- diff(range(scores$PC2, finite = TRUE))

  scores$shadow_PC1 <- scores$PC1 + (0.0045 * x_span)
  scores$shadow_PC2 <- scores$PC2 - (0.0045 * y_span)
  scores
}

mrna_pca <- compute_pca_scores(mrna_matrix, "mRNA")
mirna_pca <- compute_pca_scores(mirna_matrix, "miRNA")
pca_scores <- rbind(mrna_pca$scores, mirna_pca$scores)

make_pca_plot <- function(
  pca_result,
  panel_label,
  layer_label,
  group_variable,
  group_order,
  group_colors,
  group_fill_colors,
  base_size = 10.5,
  fixed_aspect = TRUE
) {
  scores <- pca_result$scores
  variance <- pca_result$variance
  scores$plot_group <- factor(
    scores[[group_variable]],
    levels = group_order
  )
  if (anyNA(scores$plot_group)) {
    stop(
      "PCA scores contain missing or unexpected ",
      group_variable,
      " values.",
      call. = FALSE
    )
  }
  scores <- add_point_depth_coordinates(scores)
  coordinate_system <- if (fixed_aspect) {
    ggplot2::coord_equal()
  } else {
    ggplot2::coord_cartesian()
  }
  ggplot2::ggplot(
    scores,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = plot_group,
      fill = plot_group
    )
  ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = plot_group),
      type = "norm",
      level = 0.68,
      geom = "polygon",
      alpha = 0.16,
      color = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = scores,
      ggplot2::aes(x = shadow_PC1, y = shadow_PC2),
      inherit.aes = FALSE,
      color = "#18212A",
      size = 2.85,
      alpha = 0.20,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      size = 2.45,
      alpha = 0.97,
      stroke = 0
    ) +
    ggplot2::scale_color_manual(
      values = group_colors,
      breaks = group_order,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = group_fill_colors,
      breaks = group_order,
      drop = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 4.5, alpha = 1)
      ),
      fill = "none"
    ) +
    ggplot2::labs(
      title = paste0(
        "(",
        panel_label,
        ")",
        "  Global ",
        layer_label,
        " expression structure"
      ),
      x = paste0("PC1 (", sprintf("%.1f", variance[[1]]), "%)"),
      y = paste0("PC2 (", sprintf("%.1f", variance[[2]]), "%)")
    ) +
    coordinate_system +
    publication_theme(base_size = base_size)
}

make_stratified_pca_plot <- function(
  pca_result,
  panel_label,
  layer_label,
  stratification_variable,
  stratification_order
) {
  scores <- pca_result$scores
  variance <- pca_result$variance
  scores$treatment <- factor(scores$treatment, levels = treatment_order)
  scores$stratum <- factor(
    scores[[stratification_variable]],
    levels = stratification_order
  )
  if (anyNA(scores$stratum)) {
    stop(
      "PCA scores contain missing or unexpected ",
      stratification_variable,
      " values.",
      call. = FALSE
    )
  }
  scores <- add_point_depth_coordinates(scores)

  ggplot2::ggplot(
    scores,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = treatment,
      fill = treatment
    )
  ) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = treatment),
      type = "norm",
      level = 0.68,
      geom = "polygon",
      alpha = 0.16,
      color = NA,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = scores,
      ggplot2::aes(x = shadow_PC1, y = shadow_PC2),
      inherit.aes = FALSE,
      color = "#18212A",
      size = 2.85,
      alpha = 0.20,
      stroke = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      size = 2.45,
      alpha = 0.97,
      stroke = 0
    ) +
    ggplot2::facet_wrap(~stratum, nrow = 1) +
    ggplot2::scale_color_manual(
      values = treatment_colors,
      breaks = treatment_order,
      drop = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = treatment_fill_colors,
      breaks = treatment_order,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = paste0(
        "(",
        panel_label,
        ")",
        "  Global ",
        layer_label,
        " expression structure"
      ),
      x = paste0("PC1 (", sprintf("%.1f", variance[[1]]), "%)"),
      y = paste0("PC2 (", sprintf("%.1f", variance[[2]]), "%)")
    ) +
    ggplot2::coord_equal() +
    publication_theme() +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(
        fill = "#EEF3F7",
        color = "#D0D8E0",
        linewidth = 0.45
      ),
      strip.text = ggplot2::element_text(
        face = "bold",
        color = "#263238",
        size = 10
      ),
      panel.spacing = grid::unit(0.8, "lines")
    )
}

make_experimental_timeline_plot <- function() {
  collection_points <- data.frame(
    x = c(0, 1, 4),
    y = c(3, 2, 1),
    treatment = factor(treatment_order, levels = treatment_order)
  )

  ggplot2::ggplot() +
    ggplot2::annotate(
      "segment",
      x = 0,
      xend = 1,
      y = 2,
      yend = 2,
      color = treatment_colors[["MCAO1hr"]],
      alpha = 0.86,
      linewidth = 7.0,
      lineend = "round"
    ) +
    ggplot2::annotate(
      "segment",
      x = 0,
      xend = 1,
      y = 1,
      yend = 1,
      color = treatment_colors[["MCAO3hr"]],
      alpha = 0.86,
      linewidth = 7.0,
      lineend = "round"
    ) +
    ggplot2::annotate(
      "segment",
      x = 1,
      xend = 4,
      y = 1,
      yend = 1,
      color = treatment_colors[["MCAO3hr"]],
      alpha = 0.28,
      linewidth = 7.0,
      lineend = "round"
    ) +
    ggplot2::annotate(
      "point",
      x = 1,
      y = 1,
      shape = 23,
      size = 3.3,
      fill = treatment_colors[["MCAO3hr"]],
      color = "#7D101C",
      stroke = 0.6
    ) +
    ggplot2::geom_point(
      data = collection_points,
      ggplot2::aes(x = x, y = y, color = treatment),
      size = 4.0,
      show.legend = FALSE
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 2,
      label = "60 min occlusion",
      color = "white",
      fontface = "bold",
      size = 4.0
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 1,
      label = "60 min occlusion",
      color = "white",
      fontface = "bold",
      size = 4.0
    ) +
    ggplot2::annotate(
      "text",
      x = 2.5,
      y = 1,
      label = "3 h recovery",
      color = "#7D101C",
      fontface = "bold",
      size = 4.5
    ) +
    ggplot2::annotate(
      "text",
      x = c(1, 3.72),
      y = c(2.28, 1.28),
      label = "Collection",
      color = "#28333D",
      fontface = "bold",
      size = 4.1
    ) +
    ggplot2::annotate(
      "text",
      x = 0.18,
      y = 3,
      label = "Collection",
      color = "#28333D",
      fontface = "bold",
      size = 4.1,
      hjust = 0
    ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 0.70,
      label = "Recanalization",
      color = "#7D101C",
      fontface = "bold",
      size = 3.9
    ) +
    ggplot2::scale_color_manual(values = treatment_colors) +
    ggplot2::scale_x_continuous(
      breaks = 0:4,
      limits = c(-0.12, 4.25),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(3, 2, 1),
      labels = treatment_order,
      limits = c(0.55, 3.5),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = "(A)  Experimental timeline",
      x = "Elapsed time (hours)",
      y = NULL
    ) +
    publication_theme(base_size = 14.5) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        face = "bold",
        color = "#18212A"
      ),
      plot.margin = ggplot2::margin(5, 7, 4, 7)
    )
}

make_cohort_summary_plot <- function() {
  group_centers <- c(0.17, 0.50, 0.83)
  cohort_groups <- data.frame(
    treatment = factor(treatment_order, levels = treatment_order),
    count = as.integer(observed_counts[treatment_order]),
    xmin = group_centers - 0.13,
    xmax = group_centers + 0.13,
    ymin = 0.92,
    ymax = 1.62,
    fill_color = unname(treatment_fill_colors[treatment_order]),
    border_color = unname(treatment_colors[treatment_order])
  )

  ggplot2::ggplot() +
    ggplot2::annotate(
      "rect",
      xmin = 0.23,
      xmax = 0.77,
      ymin = 2.35,
      ymax = 3.13,
      fill = "#F4F7FA",
      color = "#52616D",
      linewidth = 0.65
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 2.86,
      label = "Matched mRNA + miRNA profiles",
      color = "#18212A",
      fontface = "bold",
      size = 4.6
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 2.57,
      label = "43 animals",
      color = "#334E5C",
      fontface = "bold",
      size = 4.8
    ) +
    ggplot2::annotate(
      "segment",
      x = 0.5,
      xend = 0.5,
      y = 2.35,
      yend = 2.06,
      color = "#65747F",
      linewidth = 0.55
    ) +
    ggplot2::annotate(
      "segment",
      x = 0.17,
      xend = 0.83,
      y = 2.06,
      yend = 2.06,
      color = "#65747F",
      linewidth = 0.55
    ) +
    ggplot2::geom_segment(
      data = data.frame(x = group_centers),
      ggplot2::aes(x = x, xend = x, y = 2.06, yend = 1.67),
      color = "#65747F",
      linewidth = 0.55,
      arrow = grid::arrow(
        length = grid::unit(0.07, "inches"),
        type = "closed"
      ),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_rect(
      data = cohort_groups,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = fill_color,
        color = border_color
      ),
      linewidth = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = cohort_groups,
      ggplot2::aes(
        x = (xmin + xmax) / 2,
        y = 1.38,
        label = treatment
      ),
      color = "#18212A",
      fontface = "bold",
      size = 4.25
    ) +
    ggplot2::geom_text(
      data = cohort_groups,
      ggplot2::aes(
        x = (xmin + xmax) / 2,
        y = 1.12,
        label = paste0("n = ", count)
      ),
      color = "#334E5C",
      fontface = "bold",
      size = 4.9
    ) +
    ggplot2::annotate(
      "text",
      x = 0.5,
      y = 0.56,
      label = "24 female · 19 male     |     17 young · 26 old",
      color = "#4C5965",
      fontface = "bold",
      size = 3.9
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_color_identity() +
    ggplot2::coord_cartesian(
      xlim = c(0, 1),
      ylim = c(0.35, 3.35),
      clip = "off"
    ) +
    ggplot2::labs(title = "(B)  Paired analysis cohort") +
    ggplot2::theme_void(base_size = 14.5, base_family = "sans") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 16.0,
        color = "#18212A",
        hjust = 0
      ),
      plot.margin = ggplot2::margin(5, 7, 4, 7)
    )
}

mrna_treatment_plot <- make_pca_plot(
  mrna_pca,
  "A",
  "mRNA",
  "treatment",
  treatment_order,
  treatment_colors,
  treatment_fill_colors
)
mirna_treatment_plot <- make_pca_plot(
  mirna_pca,
  "B",
  "miRNA",
  "treatment",
  treatment_order,
  treatment_colors,
  treatment_fill_colors
)

save_composite_png(
  mrna_treatment_plot,
  mirna_treatment_plot,
  "figure1_paired_data_overview",
  12.8,
  5.4,
  component_output_dir
)

experimental_timeline_plot <- make_experimental_timeline_plot()
cohort_summary_plot <- make_cohort_summary_plot()
mrna_figure1_plot <- make_pca_plot(
  mrna_pca,
  "C",
  "mRNA",
  "treatment",
  treatment_order,
  treatment_colors,
  treatment_fill_colors,
  base_size = 14.5,
  fixed_aspect = FALSE
) +
  ggplot2::theme(plot.margin = ggplot2::margin(5, 7, 4, 7))
mirna_figure1_plot <- make_pca_plot(
  mirna_pca,
  "D",
  "miRNA",
  "treatment",
  treatment_order,
  treatment_colors,
  treatment_fill_colors,
  base_size = 14.5,
  fixed_aspect = FALSE
) +
  ggplot2::theme(plot.margin = ggplot2::margin(5, 7, 4, 7))

save_four_panel_png(
  experimental_timeline_plot,
  cohort_summary_plot,
  mrna_figure1_plot,
  mirna_figure1_plot,
  "figure1_experimental_design_and_paired_pca",
  13.2,
  8.8,
  main_output_dir
)

mrna_sex_stratified_plot <- make_stratified_pca_plot(
  mrna_pca,
  "A",
  "mRNA",
  "sex",
  sex_order
)
mirna_sex_stratified_plot <- make_stratified_pca_plot(
  mirna_pca,
  "B",
  "miRNA",
  "sex",
  sex_order
)

save_stacked_png(
  mrna_sex_stratified_plot,
  mirna_sex_stratified_plot,
  "figure_s1_pca_stratified_by_sex",
  12.8,
  9.5,
  supplementary_output_dir
)

mrna_age_stratified_plot <- make_stratified_pca_plot(
  mrna_pca,
  "A",
  "mRNA",
  "age_group",
  age_order
)
mirna_age_stratified_plot <- make_stratified_pca_plot(
  mirna_pca,
  "B",
  "miRNA",
  "age_group",
  age_order
)

save_stacked_png(
  mrna_age_stratified_plot,
  mirna_age_stratified_plot,
  "figure_s2_pca_stratified_by_age",
  12.8,
  9.5,
  supplementary_output_dir
)

write.csv(
  pca_scores,
  file.path(source_data_output_dir, "figure1_pca_scores.csv"),
  row.names = FALSE,
  na = ""
)

summary_lines <- c(
  "Manuscript Figure 1 paired-cohort PCA summary",
  "",
  paste0("Paired animals: ", nrow(manifest)),
  paste0(
    "Treatment counts: ",
    paste(
      names(observed_counts),
      as.integer(observed_counts),
      sep = "=",
      collapse = "; "
    )
  ),
  paste0(
    "Sex counts: ",
    paste(
      names(table(manifest$sex)),
      as.integer(table(manifest$sex)),
      sep = "=",
      collapse = "; "
    )
  ),
  paste0(
    "Age-group counts: ",
    paste(
      names(table(manifest$age_group)),
      as.integer(table(manifest$age_group)),
      sep = "=",
      collapse = "; "
    )
  ),
  paste0("mRNA features: ", nrow(mrna_matrix)),
  paste0("Reliably detected mature miRNA features: ", nrow(mirna_matrix)),
  paste0(
    "mRNA variance: PC1=",
    sprintf("%.2f%%", mrna_pca$variance[[1]]),
    "; PC2=",
    sprintf("%.2f%%", mrna_pca$variance[[2]])
  ),
  paste0(
    "miRNA variance: PC1=",
    sprintf("%.2f%%", mirna_pca$variance[[1]]),
    "; PC2=",
    sprintf("%.2f%%", mirna_pca$variance[[2]])
  ),
  paste0(
    "Treatment palette: Sham=",
    treatment_colors[["Sham"]],
    "; MCAO1hr=",
    treatment_colors[["MCAO1hr"]],
    "; MCAO3hr=",
    treatment_colors[["MCAO3hr"]]
  ),
  paste0(
    "Sex-stratified panels retain the treatment palette within separate ",
    "Female and Male facets."
  ),
  paste0(
    "Age-stratified panels retain the treatment palette within separate ",
    "Young and Old facets."
  ),
  "",
  paste0(
    "Interpretation boundary: PCA is descriptive and does not test ",
    "treatment effects."
  )
)
writeLines(
  summary_lines,
  file.path(source_data_output_dir, "figure1_summary.txt")
)

message(
  "Wrote publication Figure 1 outputs under: ",
  figure_root
)
