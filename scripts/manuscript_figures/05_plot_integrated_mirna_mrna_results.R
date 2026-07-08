#!/usr/bin/env Rscript

# Create the candidate-selection and network panels for manuscript Figure 4.
#
# This script visualizes the saved late-integration evidence summary. It does
# not rerun differential expression, target matching, or association testing.
# Figure 4A shows how target-supported miRNA--mRNA pairs were narrowed by
# increasingly stringent, prespecified evidence criteria. The final stage is
# intentionally shown as zero because no adjusted association met pair-level
# BH FDR < 0.10. Figure 4B shows only the ten directionally compatible but
# association-unconfirmed candidates, avoiding a large undirected network.
#
# Input:
#   results/multiomics/integrated_pair_evidence/
#     mirna_mrna_integrated_pair_evidence_summary.csv
#     mirna_mrna_directionally_compatible_unconfirmed_pairs.csv
#
# Outputs:
#   results/manuscript/figures/
#     components/figure4a_integrated_candidate_selection.png
#     components/figure4b_integrated_candidate_network.png
#     components/figure4c_inverse_pair_adjusted_associations.png
#     main/figure4_integrated_mirna_mrna.png
#     source_data/figure4a_integrated_candidate_selection.csv
#     source_data/figure4b_integrated_candidate_network.csv
#     source_data/figure4c_inverse_pair_adjusted_associations.csv
#
# Usage:
#   Rscript scripts/manuscript_figures/05_plot_integrated_mirna_mrna_results.R

required_packages <- c("ggplot2", "png")
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

input_file <- file.path(
  "results",
  "multiomics",
  "integrated_pair_evidence",
  "mirna_mrna_integrated_pair_evidence_summary.csv"
)
network_input_file <- file.path(
  "results",
  "multiomics",
  "integrated_pair_evidence",
  "mirna_mrna_directionally_compatible_unconfirmed_pairs.csv"
)
figure_root <- file.path("results", "manuscript", "figures")
main_output_dir <- file.path(figure_root, "main")
component_output_dir <- file.path(figure_root, "components")
source_data_output_dir <- file.path(figure_root, "source_data")
figure_file <- file.path(
  component_output_dir,
  "figure4a_integrated_candidate_selection.png"
)
source_data_file <- file.path(
  source_data_output_dir,
  "figure4a_integrated_candidate_selection.csv"
)
network_figure_file <- file.path(
  component_output_dir,
  "figure4b_integrated_candidate_network.png"
)
network_source_data_file <- file.path(
  source_data_output_dir,
  "figure4b_integrated_candidate_network.csv"
)
association_figure_file <- file.path(
  component_output_dir,
  "figure4c_inverse_pair_adjusted_associations.png"
)
association_source_data_file <- file.path(
  source_data_output_dir,
  "figure4c_inverse_pair_adjusted_associations.csv"
)
main_figure_file <- file.path(
  main_output_dir,
  "figure4_integrated_mirna_mrna.png"
)

missing_inputs <- c(input_file, network_input_file)[
  !file.exists(c(input_file, network_input_file))
]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing Figure 4 input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(component_output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_output_dir, recursive = TRUE, showWarnings = FALSE)

evidence <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "mirna_id",
  "gene_symbol",
  "contrast",
  "inverse_logFC_direction",
  "directionally_compatible_unconfirmed",
  "full_association_fdr_supported",
  "target_evidence_category",
  "mirna_logFC",
  "mrna_logFC",
  "adjusted_beta",
  "adjusted_SE",
  "adjusted_ci_lower",
  "adjusted_ci_upper",
  "adjusted_P.Value",
  "adjusted_adj.P.Val"
)
missing_columns <- setdiff(required_columns, names(evidence))
if (length(missing_columns) > 0L) {
  stop(
    "Integrated evidence input lacks required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

as_flag <- function(x) {
  as.character(x) %in% c("TRUE", "True", "true", "1")
}

selected_count <- nrow(evidence)
inverse_count <- sum(as_flag(evidence$inverse_logFC_direction))
compatible_count <- sum(
  as_flag(evidence$directionally_compatible_unconfirmed)
)
fdr_supported_count <- sum(
  as_flag(evidence$full_association_fdr_supported)
)
unique_mirna_count <- length(unique(evidence$mirna_id))
unique_gene_count <- length(unique(evidence$gene_symbol))
observed_contrasts <- unique(evidence$contrast)

expected_counts <- c(
  selected = 56L,
  inverse = 13L,
  compatible = 10L,
  fdr_supported = 0L,
  mirnas = 6L,
  genes = 53L
)
observed_counts <- c(
  selected = selected_count,
  inverse = inverse_count,
  compatible = compatible_count,
  fdr_supported = fdr_supported_count,
  mirnas = unique_mirna_count,
  genes = unique_gene_count
)
if (!identical(as.integer(observed_counts), as.integer(expected_counts))) {
  stop(
    "Figure 4A evidence counts differ from the validated analysis summary. ",
    "Observed: ",
    paste(names(observed_counts), observed_counts, sep = "=", collapse = ", "),
    call. = FALSE
  )
}
if (
  length(observed_contrasts) != 1L ||
    observed_contrasts != "MCAO3hr_vs_MCAO1hr"
) {
  stop(
    "Figure 4A expects all selected pairs to come from ",
    "MCAO3hr_vs_MCAO1hr.",
    call. = FALSE
  )
}

stage_data <- data.frame(
  stage_order = 1:4,
  stage = c(
    "DE-supported target pairs",
    "Inverse differential-expression direction",
    "Inverse DE plus negative adjusted association",
    "Pair-level FDR-supported associations"
  ),
  count = c(
    selected_count,
    inverse_count,
    compatible_count,
    fdr_supported_count
  ),
  retained_from_previous_percent = c(
    100,
    100 * inverse_count / selected_count,
    100 * compatible_count / inverse_count,
    0
  ),
  criterion = c(
    paste0(
      "miRNA and mRNA FDR < 0.10; target evidence; ",
      unique_mirna_count,
      " miRNAs and ",
      unique_gene_count,
      " genes"
    ),
    "Opposite miRNA and mRNA log2 fold-change directions",
    "Negative coefficient after adjustment for treatment, sex, and age",
    "Adjusted association BH FDR < 0.10"
  ),
  interpretation = c(
    "DE-supported target candidates",
    "Compatible with a simple repression model",
    "Directionally compatible but association-unconfirmed",
    "No statistically confirmed regulatory coupling"
  ),
  stringsAsFactors = FALSE
)
write.csv(stage_data, source_data_file, row.names = FALSE, na = "")

network_evidence <- read.csv(
  network_input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_network_columns <- c(
  "mirna_id",
  "gene_symbol",
  "mirna_logFC",
  "mrna_logFC",
  "target_evidence_category",
  "evidence_sources",
  "inverse_logFC_direction",
  "negative_adjusted_association",
  "fully_adjusted_all_43_beta",
  "fully_adjusted_all_43_partial_r",
  "fully_adjusted_all_43_P.Value",
  "fully_adjusted_all_43_adj.P.Val",
  "contrast"
)
missing_network_columns <- setdiff(
  required_network_columns,
  names(network_evidence)
)
if (length(missing_network_columns) > 0L) {
  stop(
    "Figure 4B input lacks required column(s): ",
    paste(missing_network_columns, collapse = ", "),
    call. = FALSE
  )
}
if (
  nrow(network_evidence) != 10L ||
    length(unique(network_evidence$mirna_id)) != 4L ||
    length(unique(network_evidence$gene_symbol)) != 9L ||
    any(!as_flag(network_evidence$inverse_logFC_direction)) ||
    any(!as_flag(network_evidence$negative_adjusted_association)) ||
    any(network_evidence$fully_adjusted_all_43_adj.P.Val < 0.10) ||
    any(network_evidence$contrast != "MCAO3hr_vs_MCAO1hr")
) {
  stop(
    "Figure 4B candidate set does not match the validated integrated summary.",
    call. = FALSE
  )
}
if (
  any(!is.finite(network_evidence$mirna_logFC)) ||
    any(!is.finite(network_evidence$mrna_logFC)) ||
    any(!is.finite(network_evidence$fully_adjusted_all_43_beta)) ||
    any(!is.finite(network_evidence$fully_adjusted_all_43_partial_r)) ||
    any(!is.finite(network_evidence$fully_adjusted_all_43_P.Value)) ||
    any(!is.finite(network_evidence$fully_adjusted_all_43_adj.P.Val))
) {
  stop("Figure 4B contains non-finite statistics.", call. = FALSE)
}

network_evidence$mirna_display <- sub(
  "^mmu-mir-",
  "miR-",
  network_evidence$mirna_id
)
network_evidence$target_evidence_display <- ifelse(
  network_evidence$target_evidence_category == "miRTarBase_only",
  "miRTarBase",
  ifelse(
    network_evidence$target_evidence_category == "TargetScanMouse_only",
    "TargetScanMouse",
    "Both resources"
  )
)
network_source_data <- network_evidence[
  ,
  c(
    "mirna_id",
    "mirna_display",
    "gene_symbol",
    "contrast",
    "mirna_logFC",
    "mrna_logFC",
    "target_evidence_display",
    "fully_adjusted_all_43_beta",
    "fully_adjusted_all_43_partial_r",
    "fully_adjusted_all_43_P.Value",
    "fully_adjusted_all_43_adj.P.Val"
  )
]
network_source_data <- network_source_data[
  order(
    network_source_data$mirna_display,
    network_source_data$gene_symbol
  ),
  ,
  drop = FALSE
]
write.csv(
  network_source_data,
  network_source_data_file,
  row.names = FALSE,
  na = ""
)

draw_stage <- function(
  y,
  width,
  height,
  fill,
  border,
  count,
  heading,
  detail,
  count_color = "white",
  text_color = "white",
  border_width = 2.5
) {
  grid::grid.roundrect(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(y, "npc"),
    width = grid::unit(width, "npc"),
    height = grid::unit(height, "npc"),
    r = grid::unit(0.025, "npc"),
    gp = grid::gpar(
      fill = fill,
      col = border,
      lwd = border_width
    )
  )
  grid::grid.text(
    label = as.character(count),
    x = grid::unit(0.5, "npc"),
    y = grid::unit(y + 0.035, "npc"),
    gp = grid::gpar(
      col = count_color,
      fontsize = 39,
      fontface = "bold",
      fontfamily = "sans"
    )
  )
  grid::grid.text(
    label = heading,
    x = grid::unit(0.5, "npc"),
    y = grid::unit(y - 0.016, "npc"),
    gp = grid::gpar(
      col = text_color,
      fontsize = 19,
      fontface = "bold",
      fontfamily = "sans"
    )
  )
  grid::grid.text(
    label = detail,
    x = grid::unit(0.5, "npc"),
    y = grid::unit(y - 0.055, "npc"),
    gp = grid::gpar(
      col = text_color,
      fontsize = 14.5,
      fontfamily = "sans"
    )
  )
}

draw_transition <- function(y_start, y_end, label) {
  grid::grid.lines(
    x = grid::unit(c(0.5, 0.5), "npc"),
    y = grid::unit(c(y_start, y_end), "npc"),
    arrow = grid::arrow(
      angle = 25,
      length = grid::unit(0.14, "inches"),
      type = "closed"
    ),
    gp = grid::gpar(
      col = "#556270",
      fill = "#556270",
      lwd = 2.6
    )
  )
  grid::grid.text(
    label = label,
    x = grid::unit(0.53, "npc"),
    y = grid::unit((y_start + y_end) / 2, "npc"),
    just = "left",
    gp = grid::gpar(
      col = "#3F4852",
      fontsize = 14,
      fontface = "bold",
      fontfamily = "sans"
    )
  )
}

grDevices::png(
  filename = figure_file,
  width = 9,
  height = 11,
  units = "in",
  res = 400,
  bg = "white"
)
grid::grid.newpage()

grid::grid.text(
  label = "(A)",
  x = grid::unit(0.035, "npc"),
  y = grid::unit(0.965, "npc"),
  just = c("left", "top"),
  gp = grid::gpar(
    col = "#111111",
    fontsize = 31,
    fontface = "bold",
    fontfamily = "sans"
  )
)
grid::grid.text(
  label = "Evidence-guided candidate prioritization",
  x = grid::unit(0.16, "npc"),
  y = grid::unit(0.959, "npc"),
  just = c("left", "top"),
  gp = grid::gpar(
    col = "#17212B",
    fontsize = 27,
    fontface = "bold",
    fontfamily = "sans"
  )
)
grid::grid.text(
  label = "Target-linked changes in MCAO3hr vs MCAO1hr",
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.902, "npc"),
  gp = grid::gpar(
    col = "#53616F",
    fontsize = 18,
    fontfamily = "sans"
  )
)

draw_stage(
  y = 0.790,
  width = 0.84,
  height = 0.155,
  fill = "#174A7E",
  border = "#174A7E",
  count = selected_count,
  heading = "DE-supported target pairs",
  detail = paste0(
    "FDR < 0.10 in both layers  |  ",
    unique_mirna_count,
    " miRNAs  |  ",
    unique_gene_count,
    " genes"
  )
)

draw_transition(
  y_start = 0.704,
  y_end = 0.651,
  label = "Inverse expression direction\n23% retained"
)

draw_stage(
  y = 0.570,
  width = 0.66,
  height = 0.145,
  fill = "#287D8E",
  border = "#287D8E",
  count = inverse_count,
  heading = "Inverse-direction pairs",
  detail = "Compatible with a simple miRNA-repression pattern"
)

draw_transition(
  y_start = 0.490,
  y_end = 0.437,
  label = "Negative adjusted association\n77% retained"
)

draw_stage(
  y = 0.355,
  width = 0.54,
  height = 0.145,
  fill = "#D3942A",
  border = "#D3942A",
  count = compatible_count,
  heading = "Direction-compatible candidates",
  detail = "Inverse DE + negative adjusted coefficient"
)

draw_transition(
  y_start = 0.274,
  y_end = 0.221,
  label = "Pair-level association FDR < 0.10"
)

draw_stage(
  y = 0.140,
  width = 0.42,
  height = 0.145,
  fill = "#FFF8F6",
  border = "#B8473B",
  count = fdr_supported_count,
  heading = "FDR-supported associations",
  detail = "No confirmed regulatory coupling",
  count_color = "#A13730",
  text_color = "#6F2924",
  border_width = 3
)

grid::grid.text(
  label = paste0(
    "Target evidence and compatible directions prioritize hypotheses;\n",
    "they do not establish direct regulation."
  ),
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.035, "npc"),
  gp = grid::gpar(
    col = "#3E4852",
    fontsize = 14.5,
    fontface = "italic",
    fontfamily = "sans"
  )
)

grDevices::dev.off()

expression_palette <- grDevices::colorRampPalette(
  c("#2B68AE", "#F7F7F7", "#C83D35")
)(201)
expression_color <- function(log2_fold_change, limit = 1) {
  clipped_value <- max(-limit, min(limit, log2_fold_change))
  palette_index <- round(
    1 + 200 * (clipped_value + limit) / (2 * limit)
  )
  expression_palette[palette_index]
}
node_text_color <- function(log2_fold_change) {
  if (abs(log2_fold_change) >= 0.55) {
    "white"
  } else {
    "#1E2832"
  }
}
format_logfc <- function(x) {
  paste0(
    "log2FC ",
    ifelse(x >= 0, "+", ""),
    sprintf("%.2f", x)
  )
}

mirna_order <- c(
  "mmu-mir-5130",
  "mmu-mir-223-3p",
  "mmu-mir-146a-5p",
  "mmu-mir-1839-3p"
)
gene_order <- c(
  "Gtf2a1",
  "Kctd3",
  "Luzp1",
  "Sestd1",
  "Cpne4",
  "Mmp16",
  "Sh3gl2",
  "Gm5464",
  "Tube1"
)
mirna_y <- stats::setNames(c(0.760, 0.535, 0.330, 0.195), mirna_order)
gene_y <- stats::setNames(
  c(0.800, 0.720, 0.640, 0.560, 0.480, 0.400, 0.320, 0.240, 0.160),
  gene_order
)

if (
  !setequal(unique(network_evidence$mirna_id), mirna_order) ||
    !setequal(unique(network_evidence$gene_symbol), gene_order)
) {
  stop(
    "Figure 4B layout identifiers differ from the validated candidate set.",
    call. = FALSE
  )
}

mirna_nodes <- unique(
  network_evidence[, c("mirna_id", "mirna_display", "mirna_logFC")]
)
mirna_nodes <- mirna_nodes[
  match(mirna_order, mirna_nodes$mirna_id),
  ,
  drop = FALSE
]
gene_nodes <- unique(
  network_evidence[, c("gene_symbol", "mrna_logFC")]
)
gene_nodes <- gene_nodes[
  match(gene_order, gene_nodes$gene_symbol),
  ,
  drop = FALSE
]

grDevices::png(
  filename = network_figure_file,
  width = 10,
  height = 11,
  units = "in",
  res = 400,
  bg = "white"
)
grid::grid.newpage()

grid::grid.text(
  label = "(B)",
  x = grid::unit(0.035, "npc"),
  y = grid::unit(0.965, "npc"),
  just = c("left", "top"),
  gp = grid::gpar(
    col = "#111111",
    fontsize = 34,
    fontface = "bold",
    fontfamily = "sans"
  )
)
grid::grid.text(
  label = "Focused miRNA-mRNA candidate network",
  x = grid::unit(0.11, "npc"),
  y = grid::unit(0.959, "npc"),
  just = c("left", "top"),
  gp = grid::gpar(
    col = "#17212B",
    fontsize = 30,
    fontface = "bold",
    fontfamily = "sans"
  )
)
grid::grid.text(
  label = "Directionally compatible candidates in MCAO3hr vs MCAO1hr",
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.905, "npc"),
  gp = grid::gpar(
    col = "#53616F",
    fontsize = 19,
    fontfamily = "sans"
  )
)

grid::grid.text(
  label = "miRNA",
  x = grid::unit(0.255, "npc"),
  y = grid::unit(0.855, "npc"),
  gp = grid::gpar(
    col = "#283640",
    fontsize = 20,
    fontface = "bold",
    fontfamily = "sans"
  )
)
grid::grid.text(
  label = "mRNA target",
  x = grid::unit(0.735, "npc"),
  y = grid::unit(0.855, "npc"),
  gp = grid::gpar(
    col = "#283640",
    fontsize = 20,
    fontface = "bold",
    fontfamily = "sans"
  )
)

# Draw target-evidence arrows first so the nodes remain visually prominent.
for (row_index in seq_len(nrow(network_evidence))) {
  edge <- network_evidence[row_index, , drop = FALSE]
  edge_source <- edge$target_evidence_display
  edge_color <- if (
    edge_source == "miRTarBase"
  ) {
    "#284F73"
  } else {
    "#775A9E"
  }
  edge_line_type <- if (
    edge_source == "miRTarBase"
  ) {
    "solid"
  } else {
    "dashed"
  }
  grid::grid.lines(
    x = grid::unit(c(0.300, 0.603), "npc"),
    y = grid::unit(
      c(
        unname(mirna_y[edge$mirna_id]),
        unname(gene_y[edge$gene_symbol])
      ),
      "npc"
    ),
    arrow = grid::arrow(
      angle = 22,
      length = grid::unit(0.11, "inches"),
      type = "closed"
    ),
    gp = grid::gpar(
      col = grDevices::adjustcolor(edge_color, alpha.f = 0.78),
      fill = grDevices::adjustcolor(edge_color, alpha.f = 0.78),
      lwd = 2.7,
      lty = edge_line_type,
      lineend = "round"
    )
  )
}

# miRNAs are circles with labels outside the node.
for (row_index in seq_len(nrow(mirna_nodes))) {
  node <- mirna_nodes[row_index, , drop = FALSE]
  y_position <- unname(mirna_y[node$mirna_id])
  grid::grid.circle(
    x = grid::unit(0.265, "npc"),
    y = grid::unit(y_position, "npc"),
    r = grid::unit(0.032, "npc"),
    gp = grid::gpar(
      fill = expression_color(node$mirna_logFC),
      col = "#FFFFFF",
      lwd = 2.4
    )
  )
  grid::grid.text(
    label = node$mirna_display,
    x = grid::unit(0.215, "npc"),
    y = grid::unit(y_position + 0.011, "npc"),
    just = "right",
    gp = grid::gpar(
      col = "#1E2832",
      fontsize = 19,
      fontface = "bold",
      fontfamily = "sans"
    )
  )
  grid::grid.text(
    label = format_logfc(node$mirna_logFC),
    x = grid::unit(0.215, "npc"),
    y = grid::unit(y_position - 0.018, "npc"),
    just = "right",
    gp = grid::gpar(
      col = "#53616F",
      fontsize = 15,
      fontfamily = "sans"
    )
  )
}

# mRNA targets are rounded rectangles with the expression estimate inside.
for (row_index in seq_len(nrow(gene_nodes))) {
  node <- gene_nodes[row_index, , drop = FALSE]
  y_position <- unname(gene_y[node$gene_symbol])
  fill_color <- expression_color(node$mrna_logFC)
  text_color <- node_text_color(node$mrna_logFC)
  grid::grid.roundrect(
    x = grid::unit(0.730, "npc"),
    y = grid::unit(y_position, "npc"),
    width = grid::unit(0.245, "npc"),
    height = grid::unit(0.064, "npc"),
    r = grid::unit(0.014, "npc"),
    gp = grid::gpar(
      fill = fill_color,
      col = "#FFFFFF",
      lwd = 2.2
    )
  )
  grid::grid.text(
    label = node$gene_symbol,
    x = grid::unit(0.730, "npc"),
    y = grid::unit(y_position + 0.010, "npc"),
    gp = grid::gpar(
      col = text_color,
      fontsize = 19,
      fontface = "bold",
      fontfamily = "sans"
    )
  )
  grid::grid.text(
    label = format_logfc(node$mrna_logFC),
    x = grid::unit(0.730, "npc"),
    y = grid::unit(y_position - 0.017, "npc"),
    gp = grid::gpar(
      col = text_color,
      fontsize = 14,
      fontfamily = "sans"
    )
  )
}

# Expression-direction legend.
grid::grid.text(
  label = "Node color: log2 fold change",
  x = grid::unit(0.245, "npc"),
  y = grid::unit(0.112, "npc"),
  gp = grid::gpar(
    col = "#283640",
    fontsize = 15.5,
    fontface = "bold",
    fontfamily = "sans"
  )
)
legend_colors <- expression_palette[seq(1, 201, length.out = 61)]
legend_x <- seq(0.105, 0.385, length.out = length(legend_colors) + 1L)
for (color_index in seq_along(legend_colors)) {
  grid::grid.rect(
    x = grid::unit(
      mean(legend_x[c(color_index, color_index + 1L)]),
      "npc"
    ),
    y = grid::unit(0.083, "npc"),
    width = grid::unit(
      diff(legend_x[c(color_index, color_index + 1L)]),
      "npc"
    ),
    height = grid::unit(0.020, "npc"),
    gp = grid::gpar(
      fill = legend_colors[color_index],
      col = NA
    )
  )
}
grid::grid.text(
  label = "lower",
  x = grid::unit(0.105, "npc"),
  y = grid::unit(0.058, "npc"),
  just = "left",
  gp = grid::gpar(
    col = "#53616F",
    fontsize = 13.5,
    fontfamily = "sans"
  )
)
grid::grid.text(
  label = "higher",
  x = grid::unit(0.385, "npc"),
  y = grid::unit(0.058, "npc"),
  just = "right",
  gp = grid::gpar(
    col = "#53616F",
    fontsize = 13.5,
    fontfamily = "sans"
  )
)

# Target-evidence legend.
grid::grid.text(
  label = "Arrow style: target evidence",
  x = grid::unit(0.705, "npc"),
  y = grid::unit(0.112, "npc"),
  gp = grid::gpar(
    col = "#283640",
    fontsize = 15.5,
    fontface = "bold",
    fontfamily = "sans"
  )
)
grid::grid.lines(
  x = grid::unit(c(0.545, 0.640), "npc"),
  y = grid::unit(c(0.083, 0.083), "npc"),
  arrow = grid::arrow(
    angle = 22,
    length = grid::unit(0.09, "inches"),
    type = "closed"
  ),
  gp = grid::gpar(
    col = "#284F73",
    fill = "#284F73",
    lwd = 2.7
  )
)
grid::grid.text(
  label = "miRTarBase",
  x = grid::unit(0.660, "npc"),
  y = grid::unit(0.083, "npc"),
  just = "left",
  gp = grid::gpar(
    col = "#384550",
    fontsize = 14,
    fontfamily = "sans"
  )
)
grid::grid.lines(
  x = grid::unit(c(0.545, 0.640), "npc"),
  y = grid::unit(c(0.055, 0.055), "npc"),
  arrow = grid::arrow(
    angle = 22,
    length = grid::unit(0.09, "inches"),
    type = "closed"
  ),
  gp = grid::gpar(
    col = "#775A9E",
    fill = "#775A9E",
    lwd = 2.7,
    lty = "dashed"
  )
)
grid::grid.text(
  label = "TargetScanMouse",
  x = grid::unit(0.660, "npc"),
  y = grid::unit(0.055, "npc"),
  just = "left",
  gp = grid::gpar(
    col = "#384550",
    fontsize = 14,
    fontfamily = "sans"
  )
)

grid::grid.text(
  label = paste0(
    "All adjusted associations were negative; none met pair-level FDR < 0.10."
  ),
  x = grid::unit(0.5, "npc"),
  y = grid::unit(0.018, "npc"),
  gp = grid::gpar(
    col = "#3E4852",
    fontsize = 14.5,
    fontface = "italic",
    fontfamily = "sans"
  )
)

grDevices::dev.off()

association_data <- evidence[
  as_flag(evidence$inverse_logFC_direction),
  ,
  drop = FALSE
]
if (
  nrow(association_data) != 13L ||
    sum(association_data$adjusted_beta < 0) != 10L ||
    sum(association_data$adjusted_beta >= 0) != 3L ||
    any(association_data$adjusted_adj.P.Val < 0.10) ||
    any(!is.finite(association_data$adjusted_beta)) ||
    any(!is.finite(association_data$adjusted_SE)) ||
    any(!is.finite(association_data$adjusted_ci_lower)) ||
    any(!is.finite(association_data$adjusted_ci_upper)) ||
    any(!is.finite(association_data$adjusted_P.Value)) ||
    any(!is.finite(association_data$adjusted_adj.P.Val)) ||
    any(association_data$adjusted_ci_lower > association_data$adjusted_beta) ||
    any(association_data$adjusted_ci_upper < association_data$adjusted_beta)
) {
  stop(
    "Figure 4C association data do not match the validated inverse-pair set.",
    call. = FALSE
  )
}

association_data$mirna_display <- sub(
  "^mmu-mir-",
  "miR-",
  association_data$mirna_id
)
association_data$pair_display <- paste(
  association_data$mirna_display,
  association_data$gene_symbol,
  sep = "  |  "
)
association_data$association_direction <- ifelse(
  association_data$adjusted_beta < 0,
  "Negative adjusted association",
  "Positive adjusted association"
)
association_data$target_evidence_display <- ifelse(
  association_data$target_evidence_category == "miRTarBase_only",
  "miRTarBase",
  ifelse(
    association_data$target_evidence_category == "TargetScanMouse_only",
    "TargetScanMouse",
    "Both resources"
  )
)
association_data$confidence_interval_crosses_zero <-
  association_data$adjusted_ci_lower <= 0 &
  association_data$adjusted_ci_upper >= 0
association_data <- association_data[
  order(association_data$adjusted_beta),
  ,
  drop = FALSE
]
association_data$y_position <- rev(seq_len(nrow(association_data)))
association_data$fdr_label <- sprintf("%.3f", association_data$adjusted_adj.P.Val)

association_source_data <- association_data[
  ,
  c(
    "mirna_id",
    "mirna_display",
    "gene_symbol",
    "pair_display",
    "contrast",
    "mirna_logFC",
    "mrna_logFC",
    "target_evidence_display",
    "adjusted_beta",
    "adjusted_SE",
    "adjusted_ci_lower",
    "adjusted_ci_upper",
    "adjusted_P.Value",
    "adjusted_adj.P.Val",
    "association_direction",
    "confidence_interval_crosses_zero"
  )
]
write.csv(
  association_source_data,
  association_source_data_file,
  row.names = FALSE,
  na = ""
)

association_plot <- ggplot2::ggplot(
  association_data,
  ggplot2::aes(y = y_position)
) +
  ggplot2::annotate(
    "rect",
    xmin = -2.8,
    xmax = 0,
    ymin = -Inf,
    ymax = Inf,
    fill = "#EFF6FB"
  ) +
  ggplot2::annotate(
    "rect",
    xmin = 0,
    xmax = 2.6,
    ymin = -Inf,
    ymax = Inf,
    fill = "#FFF7EF"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "#253746",
    linewidth = 1.1
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = adjusted_ci_lower,
      xend = adjusted_ci_upper,
      yend = y_position
    ),
    color = "#60717E",
    linewidth = 1.15,
    lineend = "round"
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = adjusted_ci_lower,
      xend = adjusted_ci_lower,
      y = y_position - 0.14,
      yend = y_position + 0.14
    ),
    color = "#60717E",
    linewidth = 0.9
  ) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = adjusted_ci_upper,
      xend = adjusted_ci_upper,
      y = y_position - 0.14,
      yend = y_position + 0.14
    ),
    color = "#60717E",
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    ggplot2::aes(
      x = adjusted_beta,
      color = association_direction,
      shape = target_evidence_display
    ),
    size = 6.2,
    stroke = 1.1
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      x = 2.34,
      label = fdr_label
    ),
    color = "#34424D",
    size = 5.8,
    fontface = "bold",
    family = "sans"
  ) +
  ggplot2::annotate(
    "text",
    x = -1.35,
    y = 13.72,
    label = "Negative coefficient",
    color = "#2A6F9B",
    size = 5.8,
    fontface = "bold",
    family = "sans"
  ) +
  ggplot2::annotate(
    "text",
    x = 0.95,
    y = 13.72,
    label = "Positive coefficient",
    color = "#B86620",
    size = 5.8,
    fontface = "bold",
    family = "sans"
  ) +
  ggplot2::annotate(
    "text",
    x = 2.34,
    y = 13.72,
    label = "FDR",
    color = "#26343E",
    size = 5.5,
    fontface = "bold",
    family = "sans"
  ) +
  ggplot2::scale_color_manual(
    values = c(
      "Negative adjusted association" = "#2A6F9B",
      "Positive adjusted association" = "#D07A27"
    ),
    name = "Adjusted direction"
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "miRTarBase" = 16,
      "TargetScanMouse" = 17,
      "Both resources" = 15
    ),
    breaks = c("miRTarBase", "TargetScanMouse"),
    name = "Target evidence"
  ) +
  ggplot2::scale_x_continuous(
    breaks = c(-2, -1, 0, 1, 2),
    limits = c(-2.8, 2.6),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_y_continuous(
    breaks = association_data$y_position,
    labels = association_data$pair_display,
    limits = c(0.45, 14.05),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "(C) Adjusted associations for inverse-direction pairs",
    subtitle = paste0(
      "mRNA expression ~ miRNA expression + treatment + sex + age"
    ),
    x = "Adjusted miRNA coefficient (95% confidence interval)",
    y = NULL,
    caption = paste0(
      "Negative estimates are directionally compatible with repression; ",
      "all confidence intervals cross zero."
    )
  ) +
  ggplot2::theme_minimal(base_family = "sans", base_size = 19) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 30,
      color = "#17212B",
      margin = ggplot2::margin(b = 7)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 20,
      color = "#53616F",
      margin = ggplot2::margin(b = 18)
    ),
    plot.caption = ggplot2::element_text(
      face = "italic",
      size = 16,
      color = "#3E4852",
      hjust = 0.5,
      margin = ggplot2::margin(t = 14)
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 17.5,
      color = "#25323C",
      margin = ggplot2::margin(r = 10)
    ),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 16.5,
      color = "#384550"
    ),
    axis.title.x = ggplot2::element_text(
      face = "bold",
      size = 20,
      color = "#26343E",
      margin = ggplot2::margin(t = 12)
    ),
    panel.grid.major.y = ggplot2::element_line(
      color = "#D9E0E5",
      linewidth = 0.55
    ),
    panel.grid.major.x = ggplot2::element_line(
      color = "#D9E0E5",
      linewidth = 0.45
    ),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 16.5,
      color = "#26343E"
    ),
    legend.text = ggplot2::element_text(
      size = 15.5,
      color = "#384550"
    ),
    legend.key.width = grid::unit(0.8, "cm"),
    plot.margin = ggplot2::margin(18, 24, 14, 18)
  )

grDevices::png(
  filename = association_figure_file,
  width = 16.5,
  height = 11,
  units = "in",
  res = 400,
  bg = "white"
)
print(association_plot)
grDevices::dev.off()

panel_a_raster <- png::readPNG(figure_file)
panel_b_raster <- png::readPNG(network_figure_file)
panel_c_raster <- png::readPNG(association_figure_file)

grDevices::png(
  filename = main_figure_file,
  width = 19,
  height = 22,
  units = "in",
  res = 400,
  bg = "white"
)
grid::grid.newpage()
figure_layout <- grid::grid.layout(
  nrow = 2,
  ncol = 2,
  widths = grid::unit(c(9, 10), "null"),
  heights = grid::unit(c(11, 11), "null")
)
grid::pushViewport(grid::viewport(layout = figure_layout))
grid::grid.raster(
  panel_a_raster,
  width = grid::unit(1, "npc"),
  height = grid::unit(1, "npc"),
  interpolate = TRUE,
  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
)
grid::grid.raster(
  panel_b_raster,
  width = grid::unit(1, "npc"),
  height = grid::unit(1, "npc"),
  interpolate = TRUE,
  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2)
)
grid::grid.raster(
  panel_c_raster,
  width = grid::unit(1, "npc"),
  height = grid::unit(1, "npc"),
  interpolate = TRUE,
  vp = grid::viewport(
    layout.pos.row = 2,
    layout.pos.col = 1:2
  )
)
grid::popViewport()
grDevices::dev.off()

message("Wrote Figure 4A: ", figure_file)
message("Wrote Figure 4A source data: ", source_data_file)
message("Wrote Figure 4B: ", network_figure_file)
message("Wrote Figure 4B source data: ", network_source_data_file)
message("Wrote Figure 4C: ", association_figure_file)
message("Wrote Figure 4C source data: ", association_source_data_file)
message("Wrote Main Figure 4: ", main_figure_file)
