#!/usr/bin/env Rscript

# Reduce redundancy and summarize the mRNA GSEA results from Script 12.
#
# This script does not rerun GSEA and does not run ORA. For each contrast and
# MSigDB collection, it uses fgsea::collapsePathways() to retain representative
# pathways and map dependent pathways to those representatives. Hallmark and
# canonical pathways remain separate because they were tested as separate FDR
# families in Script 12.
#
# The script also:
#   - expands leading-edge genes for the representative pathways;
#   - counts genes shared across representative leading edges;
#   - compares representative pathways across all three contrasts;
#   - creates a redundancy-reduced dot plot and cross-contrast NES heatmap.
#
# Inputs:
#   results/multiomics/enrichment_results/gsea/
#     combined_gsea_results.csv
#     significant_gsea_results_fdr_0.05.csv
#   results/multiomics/enrichment_inputs/gsea/
#     gsea_ranked_entrez_<contrast>.csv
#   resources/pathway_gene_sets/msigdb_2025.1.Mm/
#     mh.all.v2025.1.Mm.entrez.gmt
#     m2.cp.v2025.1.Mm.entrez.gmt
#
# Outputs:
#   results/multiomics/enrichment_results/gsea/redundancy_reduced/
#     representative_pathways.csv
#     pathway_redundancy_map.csv
#     representative_leading_edge_genes.csv
#     leading_edge_driver_gene_summary.csv
#     cross_contrast_representative_pathways.csv
#     redundancy_reduction_summary.csv
#     figures/representative_pathway_summary.png
#     figures/cross_contrast_pathway_heatmap.png
#     gsea_redundancy_summary.txt
#
# Usage:
#   Rscript scripts/multiomics/13_summarize_mrna_gsea_results.R

required_packages <- c("data.table", "fgsea", "ggplot2")
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

contrast_order <- c(
  "MCAO1hr_vs_Sham",
  "MCAO3hr_vs_Sham",
  "MCAO3hr_vs_MCAO1hr"
)
collection_order <- c("MH", "M2_CP")
gsea_result_dir <- file.path(
  "results", "multiomics", "enrichment_results", "gsea"
)
combined_result_file <- file.path(
  gsea_result_dir, "combined_gsea_results.csv"
)
significant_result_file <- file.path(
  gsea_result_dir, "significant_gsea_results_fdr_0.05.csv"
)
rank_input_dir <- file.path(
  "results", "multiomics", "enrichment_inputs", "gsea"
)
rank_files <- setNames(
  file.path(
    rank_input_dir,
    paste0("gsea_ranked_entrez_", contrast_order, ".csv")
  ),
  contrast_order
)
msigdb_dir <- file.path(
  "resources", "pathway_gene_sets", "msigdb_2025.1.Mm"
)
gmt_files <- setNames(
  file.path(
    msigdb_dir,
    c(
      "mh.all.v2025.1.Mm.entrez.gmt",
      "m2.cp.v2025.1.Mm.entrez.gmt"
    )
  ),
  collection_order
)

output_dir <- file.path(gsea_result_dir, "redundancy_reduced")
figure_dir <- file.path(output_dir, "figures")
representative_output_file <- file.path(
  output_dir, "representative_pathways.csv"
)
redundancy_map_output_file <- file.path(
  output_dir, "pathway_redundancy_map.csv"
)
leading_edge_output_file <- file.path(
  output_dir, "representative_leading_edge_genes.csv"
)
driver_output_file <- file.path(
  output_dir, "leading_edge_driver_gene_summary.csv"
)
cross_contrast_output_file <- file.path(
  output_dir, "cross_contrast_representative_pathways.csv"
)
reduction_summary_output_file <- file.path(
  output_dir, "redundancy_reduction_summary.csv"
)
representative_figure_file <- file.path(
  figure_dir, "representative_pathway_summary.png"
)
cross_contrast_figure_file <- file.path(
  figure_dir, "cross_contrast_pathway_heatmap.png"
)
run_summary_output_file <- file.path(
  output_dir, "gsea_redundancy_summary.txt"
)

fdr_threshold <- 0.05
dependency_p_value_threshold <- 0.05
dependency_permutations <- 200L
random_seed <- 20260728L

required_inputs <- c(
  combined_result_file,
  significant_result_file,
  unname(rank_files),
  unname(gmt_files)
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_gmt <- function(path) {
  fields <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  if (any(lengths(fields) < 3L)) {
    stop("Malformed GMT line(s) in: ", path, call. = FALSE)
  }
  pathway_names <- vapply(fields, `[[`, character(1), 1L)
  if (any(pathway_names == "") || anyDuplicated(pathway_names)) {
    stop(
      "GMT pathway names are missing or duplicated in: ",
      path,
      call. = FALSE
    )
  }
  pathways <- lapply(
    fields,
    function(x) unique(x[-c(1L, 2L)])
  )
  pathways <- lapply(pathways, function(x) x[x != "" & !is.na(x)])
  names(pathways) <- pathway_names
  pathways
}

read_ranking <- function(path, contrast) {
  data <- read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c(
    "ENTREZID", "SYMBOL", "rank_metric", "rank_position"
  )
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      basename(path),
      " is missing column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  data$ENTREZID <- trimws(as.character(data$ENTREZID))
  data$SYMBOL <- trimws(as.character(data$SYMBOL))
  data$rank_metric <- as.numeric(data$rank_metric)
  if (
    nrow(data) != 20082L ||
      anyNA(data$ENTREZID) ||
      any(data$ENTREZID == "") ||
      anyDuplicated(data$ENTREZID) ||
      any(!is.finite(data$rank_metric)) ||
      !identical(data$rank_position, seq_len(nrow(data))) ||
      any(diff(data$rank_metric) > 0)
  ) {
    stop(contrast, " GSEA ranking failed validation.", call. = FALSE)
  }
  data
}

resolve_exact_rank_ties <- function(statistics) {
  tie_mask <- duplicated(statistics) |
    duplicated(statistics, fromLast = TRUE)
  if (!any(tie_mask)) {
    return(statistics)
  }
  positive_gaps <- diff(sort(unique(statistics)))
  positive_gaps <- positive_gaps[positive_gaps > 0]
  tie_groups <- split(
    which(tie_mask),
    match(statistics[tie_mask], unique(statistics))
  )
  maximum_tie_size <- max(lengths(tie_groups))
  tie_step <- min(positive_gaps) /
    (maximum_tie_size + 1L) *
    1e-3
  adjusted_statistics <- statistics
  for (indices in tie_groups) {
    adjusted_statistics[indices] <- statistics[indices] -
      (seq_along(indices) - 1L) * tie_step
  }
  if (
    anyDuplicated(adjusted_statistics) ||
      any(diff(adjusted_statistics) > 0)
  ) {
    stop(
      "Deterministic GSEA tie resolution failed.",
      call. = FALSE
    )
  }
  adjusted_statistics
}

split_collapsed_values <- function(x) {
  if (is.na(x) || x == "") {
    return(character(0))
  }
  strsplit(x, ";", fixed = TRUE)[[1L]]
}

collapse_values <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != ""]
  paste(x, collapse = ";")
}

shorten_pathway_name <- function(x, width = 62L) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- sub("^BIOCARTA_", "", x)
  x <- sub("^WP_", "", x)
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- tools::toTitleCase(tolower(x))
  ifelse(
    nchar(x) > width,
    paste0(substr(x, 1L, width - 3L), "..."),
    x
  )
}

combined_results <- read.csv(
  combined_result_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
significant_results <- read.csv(
  significant_result_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_result_columns <- c(
  "contrast", "collection_id", "collection_name", "pathway",
  "measured_gene_count", "enrichment_score",
  "normalized_enrichment_score", "p_value", "FDR",
  "significant_fdr_0_05", "leading_edge_entrez_ids",
  "leading_edge_gene_symbols"
)
if (
  length(setdiff(required_result_columns, names(combined_results))) > 0L ||
    length(setdiff(
      required_result_columns,
      names(significant_results)
    )) > 0L ||
    nrow(combined_results) != 3273L ||
    nrow(significant_results) != 1474L ||
    any(significant_results$FDR >= fdr_threshold) ||
    anyDuplicated(combined_results[
      c("contrast", "collection_id", "pathway")
    ])
) {
  stop("Script 12 GSEA results failed validation.", call. = FALSE)
}

pathway_resources <- lapply(gmt_files, read_gmt)
rankings <- lapply(
  contrast_order,
  function(contrast) read_ranking(rank_files[[contrast]], contrast)
)
names(rankings) <- contrast_order

redundancy_maps <- list()
representative_tables <- list()
reduction_summaries <- list()
map_index <- 0L
representative_index <- 0L
summary_index <- 0L

for (contrast_index in seq_along(contrast_order)) {
  contrast <- contrast_order[contrast_index]
  ranking <- rankings[[contrast]]
  rank_statistics <- setNames(
    ranking$rank_metric,
    ranking$ENTREZID
  )
  rank_statistics <- resolve_exact_rank_ties(rank_statistics)

  for (collection_index in seq_along(collection_order)) {
    collection_id <- collection_order[collection_index]
    subset <- significant_results[
      significant_results$contrast == contrast &
        significant_results$collection_id == collection_id,
      ,
      drop = FALSE
    ]
    subset <- subset[
      order(
        subset$FDR,
        subset$p_value,
        -abs(subset$normalized_enrichment_score),
        subset$pathway
      ),
      ,
      drop = FALSE
    ]
    if (nrow(subset) == 0L) {
      stop(
        "No significant pathways for ",
        contrast,
        " / ",
        collection_id,
        ".",
        call. = FALSE
      )
    }
    collapse_input <- data.table::data.table(
      pathway = subset$pathway,
      pval = subset$p_value,
      padj = subset$FDR,
      ES = subset$enrichment_score,
      NES = subset$normalized_enrichment_score,
      size = subset$measured_gene_count
    )
    set.seed(
      random_seed +
        contrast_index * 100L +
        collection_index
    )
    collapsed <- fgsea::collapsePathways(
      fgseaRes = collapse_input,
      pathways = pathway_resources[[collection_id]],
      stats = rank_statistics,
      pval.threshold = dependency_p_value_threshold,
      nperm = dependency_permutations,
      gseaParam = 1
    )
    parent_pathways <- collapsed$parentPathways
    if (
      length(parent_pathways) != nrow(subset) ||
        !setequal(names(parent_pathways), subset$pathway)
    ) {
      stop(
        "Pathway dependency mapping is incomplete for ",
        contrast,
        " / ",
        collection_id,
        ".",
        call. = FALSE
      )
    }
    map <- subset[
      match(names(parent_pathways), subset$pathway),
      ,
      drop = FALSE
    ]
    resolve_terminal_parent <- function(pathway) {
      current <- pathway
      visited <- character(0)
      while (!is.na(parent_pathways[current])) {
        if (current %in% visited) {
          stop(
            "A cycle was found in the pathway dependency map.",
            call. = FALSE
          )
        }
        visited <- c(visited, current)
        current <- unname(parent_pathways[current])
      }
      current
    }
    map$representative_pathway <- vapply(
      names(parent_pathways),
      resolve_terminal_parent,
      character(1)
    )
    map$is_representative <-
      map$pathway == map$representative_pathway
    map$redundancy_status <- ifelse(
      map$is_representative,
      "representative",
      "dependent"
    )
    map <- map[
      ,
      c(
        "contrast", "collection_id", "collection_name",
        "pathway", "representative_pathway",
        "is_representative", "redundancy_status",
        "normalized_enrichment_score", "p_value", "FDR",
        "enrichment_direction", "measured_gene_count",
        "leading_edge_entrez_ids", "leading_edge_gene_symbols"
      ),
      drop = FALSE
    ]
    map_index <- map_index + 1L
    redundancy_maps[[map_index]] <- map

    representative_rows <- map[map$is_representative, , drop = FALSE]
    member_split <- split(map$pathway, map$representative_pathway)
    representative_rows$pathway_group_size <- as.integer(
      lengths(member_split)[representative_rows$pathway]
    )
    representative_rows$dependent_pathway_count <-
      representative_rows$pathway_group_size - 1L
    representative_rows$pathway_group_members <- vapply(
      representative_rows$pathway,
      function(pathway) collapse_values(member_split[[pathway]]),
      character(1)
    )
    representative_rows <- representative_rows[
      order(
        representative_rows$FDR,
        representative_rows$p_value,
        -abs(representative_rows$normalized_enrichment_score),
        representative_rows$pathway
      ),
      ,
      drop = FALSE
    ]
    row.names(representative_rows) <- NULL
    representative_index <- representative_index + 1L
    representative_tables[[representative_index]] <-
      representative_rows

    summary_index <- summary_index + 1L
    reduction_summaries[[summary_index]] <- data.frame(
      contrast = contrast,
      collection_id = collection_id,
      collection_name = subset$collection_name[1L],
      significant_pathways_before_reduction = nrow(subset),
      representative_pathways_after_reduction =
        nrow(representative_rows),
      dependent_pathways_collapsed =
        nrow(subset) - nrow(representative_rows),
      percent_reduction = round(
        100 * (
          nrow(subset) - nrow(representative_rows)
        ) / nrow(subset),
        2
      ),
      positive_representative_pathways = sum(
        representative_rows$normalized_enrichment_score > 0
      ),
      negative_representative_pathways = sum(
        representative_rows$normalized_enrichment_score < 0
      ),
      stringsAsFactors = FALSE
    )
  }
}

redundancy_map <- do.call(rbind, redundancy_maps)
row.names(redundancy_map) <- NULL
representative_pathways <- do.call(rbind, representative_tables)
row.names(representative_pathways) <- NULL
reduction_summary <- do.call(rbind, reduction_summaries)
row.names(reduction_summary) <- NULL

leading_edge_rows <- list()
leading_edge_index <- 0L
for (row_index in seq_len(nrow(representative_pathways))) {
  row <- representative_pathways[row_index, , drop = FALSE]
  gene_ids <- split_collapsed_values(row$leading_edge_entrez_ids)
  gene_symbols <- split_collapsed_values(
    row$leading_edge_gene_symbols
  )
  if (length(gene_ids) != length(gene_symbols)) {
    stop(
      "Leading-edge Entrez IDs and symbols have different lengths.",
      call. = FALSE
    )
  }
  if (length(gene_ids) == 0L) {
    next
  }
  leading_edge_index <- leading_edge_index + 1L
  leading_edge_rows[[leading_edge_index]] <- data.frame(
    contrast = row$contrast,
    collection_id = row$collection_id,
    collection_name = row$collection_name,
    representative_pathway = row$pathway,
    normalized_enrichment_score =
      row$normalized_enrichment_score,
    FDR = row$FDR,
    enrichment_direction = row$enrichment_direction,
    gene_entrez_id = gene_ids,
    gene_symbol = gene_symbols,
    stringsAsFactors = FALSE
  )
}
representative_leading_edge_genes <- do.call(
  rbind,
  leading_edge_rows
)
row.names(representative_leading_edge_genes) <- NULL

driver_split <- split(
  seq_len(nrow(representative_leading_edge_genes)),
  paste(
    representative_leading_edge_genes$contrast,
    representative_leading_edge_genes$collection_id,
    representative_leading_edge_genes$gene_entrez_id,
    sep = "|||"
  )
)
driver_summary <- do.call(
  rbind,
  lapply(
    driver_split,
    function(indices) {
      data <- representative_leading_edge_genes[
        indices,
        ,
        drop = FALSE
      ]
      data.frame(
        contrast = data$contrast[1L],
        collection_id = data$collection_id[1L],
        collection_name = data$collection_name[1L],
        gene_entrez_id = data$gene_entrez_id[1L],
        gene_symbol = data$gene_symbol[1L],
        representative_pathway_count = length(unique(
          data$representative_pathway
        )),
        representative_pathways = collapse_values(
          data$representative_pathway
        ),
        positive_nes_pathway_count = length(unique(
          data$representative_pathway[
            data$normalized_enrichment_score > 0
          ]
        )),
        negative_nes_pathway_count = length(unique(
          data$representative_pathway[
            data$normalized_enrichment_score < 0
          ]
        )),
        minimum_pathway_FDR = min(data$FDR),
        maximum_absolute_NES = max(abs(
          data$normalized_enrichment_score
        )),
        stringsAsFactors = FALSE
      )
    }
  )
)
driver_summary <- driver_summary[
  order(
    driver_summary$contrast,
    driver_summary$collection_id,
    -driver_summary$representative_pathway_count,
    driver_summary$minimum_pathway_FDR,
    driver_summary$gene_symbol
  ),
  ,
  drop = FALSE
]
row.names(driver_summary) <- NULL

representative_keys <- unique(
  representative_pathways[c("collection_id", "pathway")]
)
representative_keys$collection_name <- representative_pathways$
  collection_name[
    match(
      paste(
        representative_keys$collection_id,
        representative_keys$pathway
      ),
      paste(
        representative_pathways$collection_id,
        representative_pathways$pathway
      )
    )
  ]
cross_contrast <- representative_keys

for (contrast in contrast_order) {
  subset <- combined_results[
    combined_results$contrast == contrast,
    ,
    drop = FALSE
  ]
  subset_key <- paste(subset$collection_id, subset$pathway)
  output_key <- paste(
    cross_contrast$collection_id,
    cross_contrast$pathway
  )
  matched_index <- match(output_key, subset_key)
  if (anyNA(matched_index)) {
    stop(
      "A representative pathway is missing from contrast: ",
      contrast,
      call. = FALSE
    )
  }
  cross_contrast[[paste0("NES__", contrast)]] <-
    subset$normalized_enrichment_score[matched_index]
  cross_contrast[[paste0("FDR__", contrast)]] <-
    subset$FDR[matched_index]
  cross_contrast[[paste0("significant__", contrast)]] <-
    subset$FDR[matched_index] < fdr_threshold
  representative_key <- paste(
    representative_pathways$collection_id[
      representative_pathways$contrast == contrast
    ],
    representative_pathways$pathway[
      representative_pathways$contrast == contrast
    ]
  )
  cross_contrast[[paste0("representative__", contrast)]] <-
    output_key %in% representative_key
}

significant_columns <- paste0("significant__", contrast_order)
representative_columns <- paste0(
  "representative__", contrast_order
)
nes_columns <- paste0("NES__", contrast_order)
fdr_columns <- paste0("FDR__", contrast_order)
cross_contrast$significant_contrast_count <- rowSums(
  cross_contrast[significant_columns]
)
cross_contrast$representative_contrast_count <- rowSums(
  cross_contrast[representative_columns]
)
cross_contrast$minimum_FDR_across_contrasts <- apply(
  cross_contrast[fdr_columns],
  1L,
  min
)
cross_contrast$positive_significant_contrast_count <- vapply(
  seq_len(nrow(cross_contrast)),
  function(row_index) {
    nes <- unlist(cross_contrast[row_index, nes_columns])
    significant <- unlist(
      cross_contrast[row_index, significant_columns]
    )
    sum(nes[significant] > 0)
  },
  integer(1)
)
cross_contrast$negative_significant_contrast_count <- vapply(
  seq_len(nrow(cross_contrast)),
  function(row_index) {
    nes <- unlist(cross_contrast[row_index, nes_columns])
    significant <- unlist(
      cross_contrast[row_index, significant_columns]
    )
    sum(nes[significant] < 0)
  },
  integer(1)
)
cross_contrast$cross_contrast_pattern <- ifelse(
  cross_contrast$significant_contrast_count == 1L,
  "significant_in_one_contrast",
  ifelse(
    cross_contrast$positive_significant_contrast_count > 0L &
      cross_contrast$negative_significant_contrast_count > 0L,
    "direction_reversal",
    ifelse(
      cross_contrast$positive_significant_contrast_count > 0L,
      "same_positive_direction",
      "same_negative_direction"
    )
  )
)
cross_contrast$significant_direction_signature <- vapply(
  seq_len(nrow(cross_contrast)),
  function(row_index) {
    nes <- unlist(cross_contrast[row_index, nes_columns])
    significant <- unlist(
      cross_contrast[row_index, significant_columns]
    )
    directions <- ifelse(nes > 0, "positive", "negative")
    collapse_values(
      paste0(
        contrast_order[significant],
        ":",
        directions[significant]
      )
    )
  },
  character(1)
)
cross_contrast <- cross_contrast[
  order(
    -cross_contrast$significant_contrast_count,
    -cross_contrast$representative_contrast_count,
    cross_contrast$minimum_FDR_across_contrasts,
    cross_contrast$collection_id,
    cross_contrast$pathway
  ),
  ,
  drop = FALSE
]
row.names(cross_contrast) <- NULL

write.csv(
  representative_pathways,
  representative_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  redundancy_map,
  redundancy_map_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  representative_leading_edge_genes,
  leading_edge_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  driver_summary,
  driver_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  cross_contrast,
  cross_contrast_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  reduction_summary,
  reduction_summary_output_file,
  row.names = FALSE,
  na = ""
)

plot_split <- split(
  representative_pathways,
  interaction(
    representative_pathways$contrast,
    representative_pathways$collection_id,
    drop = TRUE
  )
)
representative_plot_data <- do.call(
  rbind,
  lapply(
    plot_split,
    function(data) {
      data <- data[
        order(
          data$FDR,
          data$p_value,
          -abs(data$normalized_enrichment_score)
        ),
        ,
        drop = FALSE
      ]
      head(data, 8L)
    }
  )
)
representative_plot_data$pathway_label <- shorten_pathway_name(
  representative_plot_data$pathway
)
representative_plot_data$plot_id <- paste(
  representative_plot_data$contrast,
  representative_plot_data$collection_id,
  representative_plot_data$pathway,
  sep = "|||"
)
representative_plot_data <- representative_plot_data[
  order(
    representative_plot_data$collection_name,
    representative_plot_data$contrast,
    representative_plot_data$normalized_enrichment_score
  ),
  ,
  drop = FALSE
]
representative_plot_data$plot_id <- factor(
  representative_plot_data$plot_id,
  levels = unique(representative_plot_data$plot_id)
)
representative_plot_labels <- setNames(
  representative_plot_data$pathway_label,
  representative_plot_data$plot_id
)
representative_plot_data$minus_log10_fdr <- -log10(
  pmax(representative_plot_data$FDR, .Machine$double.xmin)
)
representative_plot_data$direction <- ifelse(
  representative_plot_data$normalized_enrichment_score > 0,
  "Positive NES",
  "Negative NES"
)

representative_plot <- ggplot2::ggplot(
  representative_plot_data,
  ggplot2::aes(
    x = normalized_enrichment_score,
    y = plot_id,
    color = direction,
    size = minus_log10_fdr
  )
) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "grey70",
    linewidth = 0.4
  ) +
  ggplot2::geom_point(alpha = 0.85) +
  ggplot2::facet_wrap(
    ggplot2::vars(collection_name, contrast),
    ncol = 3,
    scales = "free_y"
  ) +
  ggplot2::scale_y_discrete(labels = representative_plot_labels) +
  ggplot2::scale_color_manual(
    values = c(
      "Negative NES" = "#2166AC",
      "Positive NES" = "#B2182B"
    )
  ) +
  ggplot2::labs(
    title = "Redundancy-reduced mRNA GSEA pathways",
    subtitle = paste(
      "Up to 8 representative pathways with the smallest FDR",
      "per contrast and collection"
    ),
    x = "Normalized enrichment score (NES)",
    y = NULL,
    color = "Direction",
    size = expression(-log[10](FDR))
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )
ggplot2::ggsave(
  representative_figure_file,
  representative_plot,
  width = 16,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)

heatmap_pathways <- head(cross_contrast, 36L)
heatmap_data <- combined_results[
  paste(combined_results$collection_id, combined_results$pathway) %in%
    paste(
      heatmap_pathways$collection_id,
      heatmap_pathways$pathway
    ),
  ,
  drop = FALSE
]
heatmap_data$representative <- paste(
  heatmap_data$contrast,
  heatmap_data$collection_id,
  heatmap_data$pathway
) %in% paste(
  representative_pathways$contrast,
  representative_pathways$collection_id,
  representative_pathways$pathway
)
heatmap_data$significant <- heatmap_data$FDR < fdr_threshold
heatmap_key_order <- paste(
  rev(heatmap_pathways$collection_id),
  rev(heatmap_pathways$pathway)
)
heatmap_data$pathway_key <- factor(
  paste(heatmap_data$collection_id, heatmap_data$pathway),
  levels = heatmap_key_order
)
heatmap_labels <- setNames(
  paste0(
    "[",
    rev(heatmap_pathways$collection_id),
    "] ",
    shorten_pathway_name(rev(heatmap_pathways$pathway), width = 54L)
  ),
  heatmap_key_order
)
heatmap_data$contrast <- factor(
  heatmap_data$contrast,
  levels = contrast_order
)

cross_contrast_plot <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(
    x = contrast,
    y = pathway_key,
    fill = normalized_enrichment_score,
    alpha = significant
  )
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::geom_point(
    data = heatmap_data[heatmap_data$representative, , drop = FALSE],
    shape = 8,
    color = "black",
    size = 2
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  ggplot2::scale_alpha_manual(
    values = c("FALSE" = 0.25, "TRUE" = 1),
    guide = "none"
  ) +
  ggplot2::scale_y_discrete(labels = heatmap_labels) +
  ggplot2::labs(
    title = "Cross-contrast direction of representative GSEA pathways",
    subtitle = paste(
      "Top 36 shared representatives; faded cells are not FDR",
      "significant and asterisks mark representative status"
    ),
    x = NULL,
    y = NULL,
    fill = "NES"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(
      angle = 25,
      hjust = 1
    ),
    legend.position = "bottom"
  )
ggplot2::ggsave(
  cross_contrast_figure_file,
  cross_contrast_plot,
  width = 11,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)

summary_lines <- c(
  "mRNA GSEA redundancy-reduction summary",
  "",
  paste0("Run date: ", Sys.Date()),
  paste0("Input FDR threshold: ", fdr_threshold),
  paste0("Input significant pathway rows: ", nrow(significant_results)),
  paste0(
    "Method: fgsea::collapsePathways; dependency p-value threshold = ",
    dependency_p_value_threshold,
    "; permutations = ",
    dependency_permutations
  ),
  paste0("Random seed base: ", random_seed),
  paste0(
    "FDR families remain separate by contrast and MSigDB collection."
  ),
  "",
  "Redundancy reduction by contrast and collection:"
)
for (row_index in seq_len(nrow(reduction_summary))) {
  row <- reduction_summary[row_index, , drop = FALSE]
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ",
      row$contrast,
      " / ",
      row$collection_name,
      ": ",
      row$significant_pathways_before_reduction,
      " significant -> ",
      row$representative_pathways_after_reduction,
      " representative; ",
      row$dependent_pathways_collapsed,
      " dependent pathways collapsed (",
      format(row$percent_reduction, nsmall = 2),
      "% reduction)"
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  paste0(
    "Representative pathway rows: ",
    nrow(representative_pathways)
  ),
  paste0(
    "Representative leading-edge pathway--gene rows: ",
    nrow(representative_leading_edge_genes)
  ),
  paste0("Leading-edge driver gene rows: ", nrow(driver_summary)),
  paste0(
    "Unique pathways representative in at least one contrast: ",
    nrow(cross_contrast)
  ),
  paste0(
    "Cross-contrast direction reversals: ",
    sum(cross_contrast$cross_contrast_pattern == "direction_reversal")
  ),
  "",
  paste0("Wrote representative pathways: ", representative_output_file),
  paste0("Wrote redundancy map: ", redundancy_map_output_file),
  paste0("Wrote representative leading edges: ", leading_edge_output_file),
  paste0("Wrote leading-edge driver summary: ", driver_output_file),
  paste0("Wrote cross-contrast summary: ", cross_contrast_output_file),
  paste0("Wrote reduction counts: ", reduction_summary_output_file),
  paste0("Wrote representative figure: ", representative_figure_file),
  paste0("Wrote cross-contrast figure: ", cross_contrast_figure_file),
  "",
  paste(
    "Interpretation boundary: collapsePathways reduces statistically",
    "dependent pathway results to representatives. Representative",
    "pathway names are concise labels, not proof that only that named",
    "process is active."
  ),
  paste(
    "These are mRNA GSEA summaries and do not establish direct miRNA",
    "regulation, causality, or mediation."
  )
)
writeLines(summary_lines, run_summary_output_file)

cat(
  "Reduced ",
  nrow(significant_results),
  " significant pathway rows to ",
  nrow(representative_pathways),
  " representative pathway rows.\n",
  sep = ""
)
cat("Wrote: ", run_summary_output_file, "\n", sep = "")
