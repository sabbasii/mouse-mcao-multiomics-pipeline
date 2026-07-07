#!/usr/bin/env Rscript

# Create the mRNA pathway panels for Supplementary Figure S3.
#
# This script visualizes the saved GSEA results from Script 12. It does not
# rerun differential expression, gene ranking, GSEA, or redundancy reduction.
# Twelve pathways are shown: two interpretable examples for each of six
# biological themes selected before plotting. A focused expression heatmap
# displays seven leading-edge genes per theme across the 43 paired animals.
#
# Input:
#   results/multiomics/enrichment_results/gsea/
#     combined_gsea_results.csv
#     leading_edge_genes.csv
#   results/mrna/analysis_ready/
#     expression_matrix_unique_gene_mapped_mrna.csv
#     transcript_cluster_annotation_unique_gene_mapped_mrna.csv
#   results/mrna/differential_expression/treatment_sex_age_limma/
#     analysis_samples_mrna.csv
#
# Outputs:
#   results/manuscript/figures/
#     supplementary/figure_s3_mrna_pathway_results.png
#     components/figure_s3a_mrna_gsea_nes_heatmap.png
#     components/figure_s3b_mrna_leading_edge_heatmap.png
#     source_data/figure_s3a_mrna_gsea_curated_pathways.csv
#     source_data/figure_s3b_mrna_leading_edge_expression.csv
#   results/manuscript/tables/
#     supplementary_table_s1_selected_leading_edge_genes.csv
#
# Usage:
#   Rscript scripts/manuscript_figures/03_plot_mrna_pathway_results.R

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Missing required R package: ggplot2", call. = FALSE)
}

input_file <- file.path(
  "results",
  "multiomics",
  "enrichment_results",
  "gsea",
  "combined_gsea_results.csv"
)
leading_edge_file <- file.path(
  "results",
  "multiomics",
  "enrichment_results",
  "gsea",
  "leading_edge_genes.csv"
)
expression_file <- file.path(
  "results",
  "mrna",
  "analysis_ready",
  "expression_matrix_unique_gene_mapped_mrna.csv"
)
annotation_file <- file.path(
  "results",
  "mrna",
  "analysis_ready",
  "transcript_cluster_annotation_unique_gene_mapped_mrna.csv"
)
sample_file <- file.path(
  "results",
  "mrna",
  "differential_expression",
  "treatment_sex_age_limma",
  "analysis_samples_mrna.csv"
)
figure_root <- file.path("results", "manuscript", "figures")
supplementary_output_dir <- file.path(
  figure_root,
  "supplementary"
)
component_output_dir <- file.path(figure_root, "components")
source_data_output_dir <- file.path(figure_root, "source_data")
table_output_dir <- file.path("results", "manuscript", "tables")
figure_file <- file.path(
  component_output_dir,
  "figure_s3a_mrna_gsea_nes_heatmap.png"
)
curated_output_file <- file.path(
  source_data_output_dir,
  "figure_s3a_mrna_gsea_curated_pathways.csv"
)
leading_edge_figure_file <- file.path(
  component_output_dir,
  "figure_s3b_mrna_leading_edge_heatmap.png"
)
leading_edge_output_file <- file.path(
  source_data_output_dir,
  "figure_s3b_mrna_leading_edge_expression.csv"
)
combined_figure_file <- file.path(
  supplementary_output_dir,
  "figure_s3_mrna_pathway_results.png"
)
leading_edge_table_file <- file.path(
  table_output_dir,
  "supplementary_table_s1_selected_leading_edge_genes.csv"
)

required_inputs <- c(
  input_file,
  leading_edge_file,
  expression_file,
  annotation_file,
  sample_file
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required pathway-figure input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}
dir.create(
  supplementary_output_dir,
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
dir.create(
  table_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

contrast_order <- c(
  "MCAO1hr_vs_Sham",
  "MCAO3hr_vs_Sham",
  "MCAO3hr_vs_MCAO1hr"
)
contrast_labels <- c(
  MCAO1hr_vs_Sham = "MCAO1hr\nvs Sham",
  MCAO3hr_vs_Sham = "MCAO3hr\nvs Sham",
  MCAO3hr_vs_MCAO1hr = "MCAO3hr\nvs MCAO1hr"
)

curated_pathways <- data.frame(
  theme_order = rep(seq_len(6L), each = 2L),
  pathway_order = seq_len(12L),
  theme = rep(
    c(
      "Innate inflammatory signaling",
      "Adaptive immune and antigen presentation",
      "Growth and stress signaling",
      "Mitochondrial and energy metabolism",
      "RNA processing and translation",
      "Extracellular-matrix remodeling"
    ),
    each = 2L
  ),
  pathway = c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "REACTOME_NEUTROPHIL_DEGRANULATION",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "REACTOME_CLASS_I_MHC_MEDIATED_ANTIGEN_PROCESSING_PRESENTATION",
    "HALLMARK_MTORC1_SIGNALING",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "WP_MRNA_PROCESSING",
    "WP_TRANSLATION_FACTORS",
    "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "REACTOME_ACTIVATION_OF_MATRIX_METALLOPROTEINASES"
  ),
  pathway_label = c(
    "TNFα/NF-κB signaling · Hallmark",
    "Neutrophil degranulation · Canonical",
    "Interferon-γ response · Hallmark",
    "Class I MHC antigen processing · Canonical",
    "mTORC1 signaling · Hallmark",
    "MYC targets · Hallmark",
    "Oxidative phosphorylation · Hallmark",
    "Fatty-acid metabolism · Hallmark",
    "mRNA processing · Canonical",
    "Translation factors · Canonical",
    "Extracellular-matrix organization · Canonical",
    "Matrix metalloproteinase activation · Canonical"
  ),
  expected_collection_id = c(
    "MH",
    "M2_CP",
    "MH",
    "M2_CP",
    "MH",
    "MH",
    "MH",
    "MH",
    "M2_CP",
    "M2_CP",
    "M2_CP",
    "M2_CP"
  ),
  stringsAsFactors = FALSE
)

gsea_results <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_columns <- c(
  "contrast",
  "collection_id",
  "collection_name",
  "pathway",
  "normalized_enrichment_score",
  "FDR"
)
missing_columns <- setdiff(required_columns, names(gsea_results))
if (length(missing_columns) > 0L) {
  stop(
    "Saved GSEA results lack required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}
if (
  nrow(gsea_results) != 3273L ||
    anyDuplicated(
      gsea_results[c("contrast", "collection_id", "pathway")]
    ) ||
    any(!is.finite(gsea_results$normalized_enrichment_score)) ||
    any(!is.finite(gsea_results$FDR))
) {
  stop("Saved mRNA GSEA results failed validation.", call. = FALSE)
}

curated_results <- merge(
  curated_pathways,
  gsea_results,
  by = "pathway",
  all.x = TRUE,
  sort = FALSE
)
curated_results <- curated_results[
  curated_results$contrast %in% contrast_order,
  ,
  drop = FALSE
]
expected_rows <- nrow(curated_pathways) * length(contrast_order)
if (
  nrow(curated_results) != expected_rows ||
    anyNA(curated_results$contrast) ||
    any(
      curated_results$collection_id !=
        curated_results$expected_collection_id
    ) ||
    anyDuplicated(
      curated_results[c("pathway", "contrast")]
    )
) {
  stop(
    "Curated pathways did not map uniquely to all three GSEA contrasts.",
    call. = FALSE
  )
}

curated_results <- curated_results[
  order(curated_results$pathway_order, match(
    curated_results$contrast,
    contrast_order
  )),
  ,
  drop = FALSE
]
curated_results$theme <- factor(
  curated_results$theme,
  levels = unique(curated_pathways$theme)
)
curated_results$pathway_label <- factor(
  curated_results$pathway_label,
  levels = rev(curated_pathways$pathway_label)
)
curated_results$contrast <- factor(
  curated_results$contrast,
  levels = contrast_order
)
curated_results$contrast_label <- factor(
  as.character(curated_results$contrast),
  levels = contrast_order,
  labels = unname(contrast_labels[contrast_order])
)
curated_results$significant_fdr_0_05 <-
  curated_results$FDR < 0.05
curated_results$cell_label <- paste0(
  sprintf("%+.2f", curated_results$normalized_enrichment_score),
  ifelse(curated_results$significant_fdr_0_05, " *", "")
)
curated_results$text_color <- ifelse(
  abs(curated_results$normalized_enrichment_score) >= 1.75,
  "white",
  "#17222C"
)

expected_direction <- c(
  1, -1, -1,
  1, 1, -1,
  1, -1, -1,
  1, -1, -1,
  1, -1, -1,
  1, -1, -1,
  1, -1, -1,
  1, -1, -1,
  1, -1, -1,
  1, -1, -1,
  -1, 1, 1,
  -1, 1, 1
)
observed_direction <- sign(
  curated_results$normalized_enrichment_score
)
if (!identical(observed_direction, expected_direction)) {
  stop(
    "At least one curated pathway direction differs from the documented result.",
    call. = FALSE
  )
}

write.csv(
  curated_results[
    ,
    c(
      "theme",
      "pathway_label",
      "pathway",
      "collection_id",
      "collection_name",
      "contrast",
      "normalized_enrichment_score",
      "FDR",
      "significant_fdr_0_05"
    )
  ],
  curated_output_file,
  row.names = FALSE
)

heatmap <- ggplot2::ggplot(
  curated_results,
  ggplot2::aes(
    x = contrast_label,
    y = pathway_label,
    fill = normalized_enrichment_score
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 1.15,
    width = 0.97,
    height = 0.92
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = cell_label,
      color = text_color
    ),
    size = 5.2,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(theme),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-3.7, 3.7),
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    name = "Normalized enrichment\nscore (NES)"
  ) +
  ggplot2::scale_color_identity() +
  ggplot2::labs(
    title = "(A) Cross-contrast mRNA pathway enrichment",
    subtitle = paste0(
      "Positive NES = higher-ranked genes; ",
      "negative NES = lower-ranked genes"
    ),
    x = NULL,
    y = NULL,
    caption = "* Collection-specific Benjamini–Hochberg FDR < 0.05"
  ) +
  ggplot2::theme_minimal(
    base_size = 15.5,
    base_family = "sans"
  ) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 20,
      color = "#17222C",
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = ggplot2::element_text(
      face = "bold",
      size = 13.5,
      color = "#50606D",
      margin = ggplot2::margin(b = 12)
    ),
    plot.caption = ggplot2::element_text(
      face = "bold",
      size = 12.5,
      color = "#394955",
      hjust = 0,
      margin = ggplot2::margin(t = 10)
    ),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 14,
      color = "#25313B",
      lineheight = 0.95,
      margin = ggplot2::margin(t = 8)
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 13.2,
      color = "#25313B",
      margin = ggplot2::margin(r = 8)
    ),
    panel.grid = ggplot2::element_blank(),
    panel.spacing.y = grid::unit(0.11, "lines"),
    strip.placement = "outside",
    strip.background = ggplot2::element_rect(
      fill = "#EAF3F8",
      color = "#C6D9E4",
      linewidth = 0.7
    ),
    strip.text.y.left = ggplot2::element_text(
      angle = 0,
      face = "bold",
      size = 12.6,
      color = "#24465A",
      hjust = 0.5,
      margin = ggplot2::margin(7, 9, 7, 9)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 13.2,
      color = "#25313B"
    ),
    legend.text = ggplot2::element_text(
      face = "bold",
      size = 12.5,
      color = "#34414D"
    ),
    legend.key.width = grid::unit(2.0, "cm"),
    plot.margin = ggplot2::margin(15, 18, 12, 15)
  )

grDevices::png(
  filename = figure_file,
  width = 13.2,
  height = 12.2,
  units = "in",
  res = 400,
  bg = "white"
)
print(heatmap)
grDevices::dev.off()

leading_edge_gene_selection <- data.frame(
  gene_order = seq_len(42L),
  theme = rep(unique(curated_pathways$theme), each = 7L),
  gene_symbol = c(
    "Ccrl2",
    "Zc3h12a",
    "Tnfaip3",
    "Il1b",
    "Ptgs2",
    "Dusp1",
    "Nfkbia",
    "H2-K1",
    "Tap1",
    "Tap2",
    "Psmb8",
    "Irf1",
    "Stat1",
    "Ube2l6",
    "Ddit4",
    "Hspd1",
    "Psmd14",
    "Psat1",
    "Hprt1",
    "Pgk1",
    "Eprs1",
    "Pdhb",
    "Etfa",
    "Mdh1",
    "Idh3b",
    "Aco2",
    "Vdac1",
    "Atp5mc1",
    "Eif4g3",
    "Eif4e",
    "Eif3b",
    "Eif4h",
    "Srsf3",
    "Srsf11",
    "Tra2b",
    "Vcan",
    "Mmp16",
    "Mmp17",
    "Col14a1",
    "Col10a1",
    "Jam2",
    "Colgalt2"
  ),
  stringsAsFactors = FALSE
)
if (
  anyDuplicated(leading_edge_gene_selection$gene_symbol) ||
    !identical(
      as.integer(table(leading_edge_gene_selection$theme)),
      rep(7L, 6L)
    )
) {
  stop(
    "Leading-edge gene selection must contain seven unique genes per theme.",
    call. = FALSE
  )
}

leading_edges <- read.csv(
  leading_edge_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_leading_edge_columns <- c(
  "contrast",
  "collection_id",
  "pathway",
  "NES",
  "FDR",
  "gene_entrez_id",
  "gene_symbol"
)
missing_leading_edge_columns <- setdiff(
  required_leading_edge_columns,
  names(leading_edges)
)
if (length(missing_leading_edge_columns) > 0L) {
  stop(
    "Saved leading-edge results lack required column(s): ",
    paste(missing_leading_edge_columns, collapse = ", "),
    call. = FALSE
  )
}
curated_leading_edges <- merge(
  leading_edges,
  curated_pathways[c("theme", "pathway")],
  by = "pathway",
  all = FALSE,
  sort = FALSE
)
correct_theme_membership <- merge(
  leading_edge_gene_selection,
  curated_leading_edges,
  by = c("theme", "gene_symbol"),
  all.x = TRUE,
  sort = FALSE
)
genes_without_membership <- unique(
  correct_theme_membership$gene_symbol[
    is.na(correct_theme_membership$pathway)
  ]
)
if (length(genes_without_membership) > 0L) {
  stop(
    "Selected gene(s) lack leading-edge membership in their theme: ",
    paste(genes_without_membership, collapse = ", "),
    call. = FALSE
  )
}

leading_edge_membership_summary <- aggregate(
  cbind(
    selected_pathway_count = correct_theme_membership$pathway,
    contrast_count = correct_theme_membership$contrast
  ) ~ theme + gene_symbol,
  data = correct_theme_membership,
  FUN = function(x) length(unique(x))
)
leading_edge_minimum_fdr <- aggregate(
  FDR ~ theme + gene_symbol,
  data = correct_theme_membership,
  FUN = min
)
names(leading_edge_minimum_fdr)[
  names(leading_edge_minimum_fdr) == "FDR"
] <- "minimum_pathway_FDR"
leading_edge_gene_selection <- merge(
  leading_edge_gene_selection,
  leading_edge_membership_summary,
  by = c("theme", "gene_symbol"),
  all.x = TRUE,
  sort = FALSE
)
leading_edge_gene_selection <- merge(
  leading_edge_gene_selection,
  leading_edge_minimum_fdr,
  by = c("theme", "gene_symbol"),
  all.x = TRUE,
  sort = FALSE
)
leading_edge_gene_selection <- leading_edge_gene_selection[
  order(leading_edge_gene_selection$gene_order),
  ,
  drop = FALSE
]
leading_edge_gene_selection$theme_order <- match(
  leading_edge_gene_selection$theme,
  unique(curated_pathways$theme)
)
leading_edge_gene_selection$gene_order_within_theme <- ave(
  leading_edge_gene_selection$gene_order,
  leading_edge_gene_selection$theme,
  FUN = seq_along
)

write.csv(
  leading_edge_gene_selection[
    ,
    c(
      "theme_order",
      "theme",
      "gene_order_within_theme",
      "gene_symbol",
      "selected_pathway_count",
      "contrast_count",
      "minimum_pathway_FDR"
    )
  ],
  leading_edge_table_file,
  row.names = FALSE
)

analysis_samples <- read.csv(
  sample_file,
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
if (
  length(missing_sample_columns) > 0L ||
    nrow(analysis_samples) != 43L ||
    anyDuplicated(analysis_samples$file_name) ||
    anyDuplicated(analysis_samples$animal_id)
) {
  stop(
    "Saved paired mRNA sample metadata failed validation.",
    call. = FALSE
  )
}
treatment_order <- c("Sham", "MCAO1hr", "MCAO3hr")
sex_order <- c("Female", "Male")
age_order <- c("Young", "Old")
analysis_samples$treatment <- factor(
  analysis_samples$treatment,
  levels = treatment_order
)
analysis_samples$sex <- factor(
  analysis_samples$sex,
  levels = sex_order
)
analysis_samples$age_group <- factor(
  analysis_samples$age_group,
  levels = age_order
)
if (
  anyNA(analysis_samples$treatment) ||
    anyNA(analysis_samples$sex) ||
    anyNA(analysis_samples$age_group)
) {
  stop(
    "Treatment, sex, or age annotation is missing for a paired animal.",
    call. = FALSE
  )
}
sample_order <- order(
  analysis_samples$treatment,
  analysis_samples$sex,
  analysis_samples$age_group,
  analysis_samples$animal_id
)
analysis_samples <- analysis_samples[
  sample_order,
  ,
  drop = FALSE
]
analysis_samples$sample_index <- seq_len(nrow(analysis_samples))

# Give each treatment the same compact visual width while retaining every
# animal. Samples remain ordered within treatment by sex, age, and animal ID.
heatmap_group_width <- 11
analysis_samples$treatment_index <- match(
  as.character(analysis_samples$treatment),
  treatment_order
)
analysis_samples$within_treatment_index <- ave(
  analysis_samples$sample_index,
  analysis_samples$treatment,
  FUN = seq_along
)
analysis_samples$treatment_sample_count <- as.integer(
  table(analysis_samples$treatment)[analysis_samples$treatment]
)
analysis_samples$heatmap_tile_width <-
  heatmap_group_width / analysis_samples$treatment_sample_count
analysis_samples$heatmap_x <-
  (analysis_samples$treatment_index - 1) * heatmap_group_width +
  (analysis_samples$within_treatment_index - 0.5) *
  analysis_samples$heatmap_tile_width

annotation <- read.csv(
  annotation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (
  !all(
    c("transcript_cluster_id", "SYMBOL") %in% names(annotation)
  ) ||
    anyDuplicated(annotation$transcript_cluster_id)
) {
  stop("Analysis-ready mRNA annotation failed validation.", call. = FALSE)
}
annotation_rows <- match(
  leading_edge_gene_selection$gene_symbol,
  annotation$SYMBOL
)
if (
  anyNA(annotation_rows) ||
    anyDuplicated(
      annotation$SYMBOL[
        annotation$SYMBOL %in%
          leading_edge_gene_selection$gene_symbol
      ]
    )
) {
  stop(
    "Selected leading-edge genes do not map uniquely to mRNA features.",
    call. = FALSE
  )
}
leading_edge_gene_selection$transcript_cluster_id <-
  annotation$transcript_cluster_id[annotation_rows]

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
    "Analysis-ready mRNA expression lacks required features or samples.",
    call. = FALSE
  )
}
expression_rows <- match(
  leading_edge_gene_selection$transcript_cluster_id,
  expression_table$transcript_cluster_id
)
if (anyNA(expression_rows)) {
  stop(
    "At least one selected leading-edge feature is absent from expression.",
    call. = FALSE
  )
}
leading_edge_expression <- as.matrix(
  expression_table[
    expression_rows,
    analysis_samples$file_name,
    drop = FALSE
  ]
)
storage.mode(leading_edge_expression) <- "numeric"
rownames(leading_edge_expression) <-
  leading_edge_gene_selection$gene_symbol
if (any(!is.finite(leading_edge_expression))) {
  stop(
    "Leading-edge expression contains non-finite values.",
    call. = FALSE
  )
}
gene_standard_deviations <- apply(
  leading_edge_expression,
  1L,
  stats::sd
)
if (any(!is.finite(gene_standard_deviations)) ||
    any(gene_standard_deviations == 0)) {
  stop(
    "At least one selected leading-edge gene has zero expression variance.",
    call. = FALSE
  )
}
scaled_expression <- t(
  scale(t(leading_edge_expression))
)

leading_edge_plot_data <- data.frame(
  gene_symbol = rep(
    rownames(leading_edge_expression),
    each = ncol(leading_edge_expression)
  ),
  sample_index = rep(
    analysis_samples$sample_index,
    times = nrow(leading_edge_expression)
  ),
  heatmap_x = rep(
    analysis_samples$heatmap_x,
    times = nrow(leading_edge_expression)
  ),
  heatmap_tile_width = rep(
    analysis_samples$heatmap_tile_width,
    times = nrow(leading_edge_expression)
  ),
  animal_id = rep(
    analysis_samples$animal_id,
    times = nrow(leading_edge_expression)
  ),
  treatment = rep(
    as.character(analysis_samples$treatment),
    times = nrow(leading_edge_expression)
  ),
  sex = rep(
    as.character(analysis_samples$sex),
    times = nrow(leading_edge_expression)
  ),
  age_group = rep(
    as.character(analysis_samples$age_group),
    times = nrow(leading_edge_expression)
  ),
  expression = as.vector(t(leading_edge_expression)),
  z_score = as.vector(t(scaled_expression)),
  stringsAsFactors = FALSE
)
gene_metadata_index <- match(
  leading_edge_plot_data$gene_symbol,
  leading_edge_gene_selection$gene_symbol
)
leading_edge_plot_data$theme <-
  leading_edge_gene_selection$theme[gene_metadata_index]
leading_edge_plot_data$gene_order <-
  leading_edge_gene_selection$gene_order[gene_metadata_index]
leading_edge_plot_data$gene_row <-
  nrow(leading_edge_gene_selection) -
  leading_edge_plot_data$gene_order + 1L
leading_edge_plot_data$display_z <- pmax(
  -2,
  pmin(2, leading_edge_plot_data$z_score)
)

write.csv(
  merge(
    leading_edge_plot_data,
    leading_edge_gene_selection[
      ,
      c(
        "gene_symbol",
        "transcript_cluster_id",
        "selected_pathway_count",
        "contrast_count",
        "minimum_pathway_FDR"
      )
    ],
    by = "gene_symbol",
    all.x = TRUE,
    sort = FALSE
  ),
  leading_edge_output_file,
  row.names = FALSE
)

annotation_colors <- c(
  Sham = "#13A62A",
  MCAO1hr = "#1749D1",
  MCAO3hr = "#E32636",
  Female = "#B45ACB",
  Male = "#3B82B5",
  Young = "#F0A928",
  Old = "#667985"
)
annotation_rows_data <- rbind(
  data.frame(
    sample_index = analysis_samples$sample_index,
    heatmap_x = analysis_samples$heatmap_x,
    heatmap_tile_width = analysis_samples$heatmap_tile_width,
    annotation_row = 48,
    category = as.character(analysis_samples$treatment)
  ),
  data.frame(
    sample_index = analysis_samples$sample_index,
    heatmap_x = analysis_samples$heatmap_x,
    heatmap_tile_width = analysis_samples$heatmap_tile_width,
    annotation_row = 46,
    category = as.character(analysis_samples$sex)
  ),
  data.frame(
    sample_index = analysis_samples$sample_index,
    heatmap_x = analysis_samples$heatmap_x,
    heatmap_tile_width = analysis_samples$heatmap_tile_width,
    annotation_row = 44,
    category = as.character(analysis_samples$age_group)
  )
)
annotation_rows_data$color <- unname(
  annotation_colors[annotation_rows_data$category]
)
if (anyNA(annotation_rows_data$color)) {
  stop("An annotation color is undefined.", call. = FALSE)
}

theme_labels <- c(
  "Innate inflammation",
  "Adaptive immunity",
  "Growth / stress",
  "Mitochondrial energy",
  "RNA processing",
  "ECM remodeling"
)
theme_rectangles <- do.call(
  rbind,
  lapply(
    seq_along(unique(curated_pathways$theme)),
    function(index) {
      theme_name <- unique(curated_pathways$theme)[index]
      gene_rows <- unique(
        leading_edge_plot_data$gene_row[
          leading_edge_plot_data$theme == theme_name
        ]
      )
      data.frame(
        theme = theme_name,
        theme_label = theme_labels[index],
        ymin = min(gene_rows) - 0.48,
        ymax = max(gene_rows) + 0.48,
        y = mean(range(gene_rows)),
        fill = if (index %% 2L == 1L) "#EAF3F8" else "#F3F7FA",
        stringsAsFactors = FALSE
      )
    }
  )
)
gene_axis_rows <- nrow(leading_edge_gene_selection) -
  leading_edge_gene_selection$gene_order + 1L
y_breaks <- c(48, 46, 44)
y_labels <- c("Treatment", "Sex", "Age")
annotation_label_data <- data.frame(
  x = rep(-0.45, length(y_breaks)),
  annotation_row = y_breaks,
  label = y_labels,
  stringsAsFactors = FALSE
)
treatment_counts <- table(analysis_samples$treatment)
treatment_centers <- (
  seq_along(treatment_order) - 0.5
) * heatmap_group_width
treatment_boundaries <- seq_len(
  length(treatment_order) - 1L
) * heatmap_group_width

leading_edge_heatmap <- ggplot2::ggplot() +
  ggplot2::geom_rect(
    data = theme_rectangles,
    ggplot2::aes(
      xmin = -18.50,
      xmax = -0.25,
      ymin = ymin,
      ymax = ymax
    ),
    fill = theme_rectangles$fill,
    color = "#CADCE6",
    linewidth = 0.55,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = theme_rectangles,
    ggplot2::aes(
      x = -9.35,
      y = y,
      label = theme_label
    ),
    color = "#24465A",
    fontface = "bold",
    size = 5.1,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_tile(
    data = leading_edge_plot_data,
    ggplot2::aes(
      x = heatmap_x,
      y = gene_row,
      fill = display_z,
      width = heatmap_tile_width
    ),
    height = 0.94,
    color = NA
  ) +
  ggplot2::geom_vline(
    xintercept = treatment_boundaries,
    color = "white",
    linewidth = 1.35
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#D73027",
    midpoint = 0,
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    name = "Gene-wise expression\nz score"
  )

for (category_name in names(annotation_colors)) {
  category_data <- annotation_rows_data[
    annotation_rows_data$category == category_name,
    ,
    drop = FALSE
  ]
  leading_edge_heatmap <- leading_edge_heatmap +
    ggplot2::geom_tile(
      data = category_data,
      ggplot2::aes(
        x = heatmap_x,
        y = annotation_row,
        width = heatmap_tile_width
      ),
      fill = annotation_colors[[category_name]],
      height = 1.30,
      color = NA,
      inherit.aes = FALSE
    )
}

leading_edge_heatmap <- leading_edge_heatmap +
  ggplot2::geom_text(
    data = annotation_label_data,
    ggplot2::aes(
      x = x,
      y = annotation_row,
      label = label
    ),
    hjust = 1,
    color = "#25313B",
    fontface = "bold",
    size = 5.2,
    inherit.aes = FALSE
  )

annotation_legend_data <- data.frame(
  x = rep(45, length(annotation_colors)),
  y = rep(0, length(annotation_colors)),
  category = factor(
    names(annotation_colors),
    levels = names(annotation_colors)
  )
)
leading_edge_heatmap <- leading_edge_heatmap +
  ggplot2::geom_point(
    data = annotation_legend_data,
    ggplot2::aes(x = x, y = y, color = category),
    shape = 15,
    size = 5.5,
    inherit.aes = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = annotation_colors,
    breaks = names(annotation_colors),
    name = "Sample annotation"
  ) +
  ggplot2::scale_x_continuous(
    breaks = treatment_centers,
    labels = names(treatment_counts),
    expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  ggplot2::coord_cartesian(
    xlim = c(-18.75, 33.05),
    ylim = c(0.45, 48.75),
    clip = "off"
  ) +
  ggplot2::labs(
    title = "(B) Leading-edge genes across paired animals",
    subtitle = paste0(
      "Expression standardized within each gene; samples ordered by\n",
      "treatment, sex, age, and animal ID"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Genes and row order: Supplementary Table S1. ",
      "Bars: treatment, sex, and age."
    )
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      order = 1,
      title.position = "top",
      barwidth = grid::unit(4.0, "cm"),
      barheight = grid::unit(0.55, "cm")
    ),
    color = ggplot2::guide_legend(
      order = 2,
      title.position = "top",
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 6)
    )
  ) +
  ggplot2::theme_minimal(
    base_size = 16,
    base_family = "sans"
  ) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 24,
      color = "#17222C",
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = ggplot2::element_text(
      face = "bold",
      size = 14,
      color = "#50606D",
      margin = ggplot2::margin(b = 10)
    ),
    plot.caption = ggplot2::element_text(
      face = "bold",
      size = 12,
      color = "#394955",
      hjust = 0,
      margin = ggplot2::margin(t = 9)
    ),
    axis.text.x = ggplot2::element_text(
      face = "bold",
      size = 17,
      color = "#25313B",
      margin = ggplot2::margin(t = 7)
    ),
    axis.text.y = ggplot2::element_text(
      face = "bold",
      size = 13.5,
      color = "#25313B",
      margin = ggplot2::margin(r = 4)
    ),
    panel.grid = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = ggplot2::element_text(
      face = "bold",
      size = 13,
      color = "#25313B"
    ),
    legend.text = ggplot2::element_text(
      face = "bold",
      size = 12,
      color = "#34414D"
    ),
    plot.margin = ggplot2::margin(8, 12, 8, 10)
  )

grDevices::png(
  filename = leading_edge_figure_file,
  width = 10.5,
  height = 7.6,
  units = "in",
  res = 400,
  bg = "white"
)
print(leading_edge_heatmap)
grDevices::dev.off()

draw_combined_pathway_figure <- function() {
  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 2L,
    ncol = 1L,
    heights = grid::unit(c(12.2, 7.6), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(
    heatmap,
    vp = grid::viewport(
      layout.pos.row = 1L,
      layout.pos.col = 1L
    )
  )
  print(
    leading_edge_heatmap,
    vp = grid::viewport(
      layout.pos.row = 2L,
      layout.pos.col = 1L,
      width = 0.86
    )
  )
  grid::popViewport()
}

grDevices::png(
  filename = combined_figure_file,
  width = 13.2,
  height = 19.8,
  units = "in",
  res = 400,
  bg = "white"
)
draw_combined_pathway_figure()
grDevices::dev.off()

message("Wrote Supplementary Figure S3A: ", figure_file)
message("Wrote Supplementary Figure S3B: ", leading_edge_figure_file)
message("Wrote combined Supplementary Figure S3: ", combined_figure_file)
message("Wrote curated GSEA values: ", curated_output_file)
message("Wrote leading-edge expression: ", leading_edge_output_file)
message("Wrote Supplementary Table S1: ", leading_edge_table_file)
