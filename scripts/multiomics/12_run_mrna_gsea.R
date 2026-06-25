#!/usr/bin/env Rscript

# Run gene set enrichment analysis (GSEA) on the complete mRNA rankings.
#
# This script performs mRNA GSEA only. It does not run over-representation
# analysis (ORA), test the 53-gene selected list, or interpret pathways as
# evidence of direct miRNA regulation.
#
# Each treatment contrast uses all 20,082 unique measured Entrez genes ranked
# by the limma moderated t-statistic prepared by Script 11. A positive
# normalized enrichment score (NES) means that a pathway is concentrated
# toward genes with positive t-statistics; a negative NES means that it is
# concentrated toward genes with negative t-statistics.
#
# Mouse MSigDB 2025.1 collections are tested separately:
#   - MH: 50 Hallmark gene sets for a compact biological overview;
#   - M2:CP: curated canonical pathways for detailed pathway coverage.
#
# The Hallmark and canonical-pathway results receive separate Benjamini-
# Hochberg FDR corrections within each contrast. Gene sets are tested after
# intersection with the measured genes, using sizes from 15 through 500.
#
# Inputs:
#   results/multiomics/enrichment_inputs/gsea/
#     gsea_ranked_entrez_<contrast>.csv
#   resources/pathway_gene_sets/msigdb_2025.1.Mm/
#     mh.all.v2025.1.Mm.entrez.gmt
#     m2.cp.v2025.1.Mm.entrez.gmt
#
# Outputs:
#   results/multiomics/enrichment_results/gsea/
#     contrast_results/<contrast>_gsea_results.csv
#     combined_gsea_results.csv
#     significant_gsea_results_fdr_0.05.csv
#     leading_edge_genes.csv
#     pathway_collection_summary.csv
#     figures/gsea_pathway_summary.png
#     gsea_run_summary.txt
#
# Usage:
#   Rscript scripts/multiomics/12_run_mrna_gsea.R

required_packages <- c("BiocParallel", "fgsea", "ggplot2")
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

msigdb_version <- "2025.1.Mm"
msigdb_dir <- file.path(
  "resources", "pathway_gene_sets", "msigdb_2025.1.Mm"
)
collection_definitions <- data.frame(
  collection_id = c("MH", "M2_CP"),
  collection_name = c(
    "Mouse Hallmark",
    "Mouse canonical pathways"
  ),
  gmt_file = file.path(
    msigdb_dir,
    c(
      "mh.all.v2025.1.Mm.entrez.gmt",
      "m2.cp.v2025.1.Mm.entrez.gmt"
    )
  ),
  expected_sha256 = c(
    paste0(
      "d66c8b1ac89e533b581651b98cfa47c660f49481564e23",
      "fe958781b744a1e60d"
    ),
    paste0(
      "879fd7e2c9869bd94aab08a04a5e0caf085baab818ed7",
      "c05017a1f690d938c80"
    )
  ),
  stringsAsFactors = FALSE
)

output_dir <- file.path(
  "results", "multiomics", "enrichment_results", "gsea"
)
contrast_output_dir <- file.path(output_dir, "contrast_results")
figure_output_dir <- file.path(output_dir, "figures")
combined_output_file <- file.path(
  output_dir, "combined_gsea_results.csv"
)
significant_output_file <- file.path(
  output_dir, "significant_gsea_results_fdr_0.05.csv"
)
leading_edge_output_file <- file.path(
  output_dir, "leading_edge_genes.csv"
)
collection_summary_file <- file.path(
  output_dir, "pathway_collection_summary.csv"
)
summary_output_file <- file.path(
  output_dir, "gsea_run_summary.txt"
)
summary_figure_file <- file.path(
  figure_output_dir, "gsea_pathway_summary.png"
)

minimum_gene_set_size <- 15L
maximum_gene_set_size <- 500L
fdr_threshold <- 0.05
random_seed <- 20260728L

required_inputs <- c(unname(rank_files), collection_definitions$gmt_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

dir.create(
  contrast_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  figure_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

sha256_file <- function(path) {
  output <- system2(
    "shasum",
    c("-a", "256", path),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Could not calculate SHA-256 for: ", path, call. = FALSE)
  }
  strsplit(output[1L], "[[:space:]]+")[[1L]][1L]
}

observed_sha256 <- vapply(
  collection_definitions$gmt_file,
  sha256_file,
  character(1)
)
if (!identical(
  unname(observed_sha256),
  collection_definitions$expected_sha256
)) {
  stop(
    "At least one Mouse MSigDB GMT checksum differs from the ",
    "documented source file.",
    call. = FALSE
  )
}
collection_definitions$sha256 <- observed_sha256

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  malformed <- lengths(fields) < 3L
  if (any(malformed)) {
    stop("Malformed GMT line(s) in: ", path, call. = FALSE)
  }
  pathway_names <- vapply(fields, `[[`, character(1), 1L)
  descriptions <- vapply(fields, `[[`, character(1), 2L)
  if (
    anyNA(pathway_names) ||
      any(pathway_names == "") ||
      anyDuplicated(pathway_names)
  ) {
    stop(
      "GMT pathway names are missing or duplicated in: ",
      path,
      call. = FALSE
    )
  }
  pathways <- lapply(fields, function(x) unique(x[-c(1L, 2L)]))
  pathways <- lapply(
    pathways,
    function(x) x[!is.na(x) & x != ""]
  )
  names(pathways) <- pathway_names
  metadata <- data.frame(
    pathway = pathway_names,
    pathway_description = descriptions,
    original_gene_count = lengths(pathways),
    stringsAsFactors = FALSE
  )
  list(pathways = pathways, metadata = metadata)
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
    stop(
      contrast,
      " GSEA ranking failed validation.",
      call. = FALSE
    )
  }
  data
}

resolve_exact_rank_ties <- function(statistics) {
  tie_mask <- duplicated(statistics) |
    duplicated(statistics, fromLast = TRUE)
  if (!any(tie_mask)) {
    return(list(
      statistics = statistics,
      tie_groups = 0L,
      tied_genes = 0L,
      maximum_adjustment = 0
    ))
  }
  unique_statistics <- sort(unique(statistics))
  positive_gaps <- diff(unique_statistics)
  positive_gaps <- positive_gaps[positive_gaps > 0]
  if (length(positive_gaps) == 0L) {
    stop(
      "All GSEA rank statistics are identical.",
      call. = FALSE
    )
  }
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
  list(
    statistics = adjusted_statistics,
    tie_groups = length(tie_groups),
    tied_genes = sum(lengths(tie_groups)),
    maximum_adjustment =
      (maximum_tie_size - 1L) * tie_step
  )
}

collapse_values <- function(x) {
  x <- unique(as.character(x))
  x <- x[!is.na(x) & x != ""]
  paste(x, collapse = ";")
}

make_empty_leading_edge <- function() {
  data.frame(
    contrast = character(0),
    collection_id = character(0),
    collection_name = character(0),
    pathway = character(0),
    NES = numeric(0),
    p_value = numeric(0),
    FDR = numeric(0),
    gene_entrez_id = character(0),
    gene_symbol = character(0),
    stringsAsFactors = FALSE
  )
}

collection_resources <- lapply(
  collection_definitions$gmt_file,
  read_gmt
)
names(collection_resources) <- collection_definitions$collection_id

all_results <- list()
all_leading_edge <- list()
collection_summaries <- list()
tie_summaries <- list()
result_index <- 0L
leading_edge_index <- 0L
summary_index <- 0L

for (contrast_index in seq_along(contrast_order)) {
  contrast <- contrast_order[contrast_index]
  ranking <- read_ranking(rank_files[[contrast]], contrast)
  rank_statistics <- setNames(
    ranking$rank_metric,
    ranking$ENTREZID
  )
  tie_resolution <- resolve_exact_rank_ties(rank_statistics)
  rank_statistics <- tie_resolution$statistics
  tie_summaries[[contrast]] <- data.frame(
    contrast = contrast,
    exact_tie_groups = tie_resolution$tie_groups,
    genes_in_exact_ties = tie_resolution$tied_genes,
    maximum_absolute_tie_adjustment =
      tie_resolution$maximum_adjustment,
    stringsAsFactors = FALSE
  )
  symbol_by_entrez <- setNames(ranking$SYMBOL, ranking$ENTREZID)
  contrast_results <- list()

  for (
    collection_index in seq_len(nrow(collection_definitions))
  ) {
    collection_id <- collection_definitions$collection_id[
      collection_index
    ]
    collection_name <- collection_definitions$collection_name[
      collection_index
    ]
    resource <- collection_resources[[collection_id]]
    set.seed(
      random_seed +
        (contrast_index * 100L) +
        collection_index
    )
    fgsea_result <- fgsea::fgseaMultilevel(
      pathways = resource$pathways,
      stats = rank_statistics,
      minSize = minimum_gene_set_size,
      maxSize = maximum_gene_set_size,
      eps = 0,
      scoreType = "std",
      BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
    )
    result <- as.data.frame(fgsea_result)
    if (nrow(result) == 0L) {
      stop(
        "No pathways were tested for ",
        contrast,
        " / ",
        collection_name,
        ".",
        call. = FALSE
      )
    }
    pathway_metadata <- resource$metadata[
      match(result$pathway, resource$metadata$pathway),
      ,
      drop = FALSE
    ]
    if (anyNA(pathway_metadata$pathway)) {
      stop(
        "Could not match GSEA results to pathway metadata.",
        call. = FALSE
      )
    }
    leading_edge_entrez <- vapply(
      result$leadingEdge,
      collapse_values,
      character(1)
    )
    leading_edge_symbols <- vapply(
      result$leadingEdge,
      function(ids) {
        collapse_values(unname(symbol_by_entrez[ids]))
      },
      character(1)
    )
    result_output <- data.frame(
      contrast = contrast,
      collection_id = collection_id,
      collection_name = collection_name,
      msigdb_version = msigdb_version,
      pathway = result$pathway,
      pathway_description =
        pathway_metadata$pathway_description,
      original_gene_count =
        pathway_metadata$original_gene_count,
      measured_gene_count = result$size,
      enrichment_score = result$ES,
      normalized_enrichment_score = result$NES,
      p_value = result$pval,
      FDR = result$padj,
      log2_error = result$log2err,
      enrichment_direction = ifelse(
        result$NES > 0,
        "positive_t_statistics",
        "negative_t_statistics"
      ),
      significant_fdr_0_05 =
        !is.na(result$padj) & result$padj < fdr_threshold,
      leading_edge_entrez_ids = leading_edge_entrez,
      leading_edge_gene_symbols = leading_edge_symbols,
      stringsAsFactors = FALSE
    )
    result_output <- result_output[
      order(
        result_output$FDR,
        result_output$p_value,
        -abs(result_output$normalized_enrichment_score),
        result_output$pathway,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]
    row.names(result_output) <- NULL

    result_index <- result_index + 1L
    all_results[[result_index]] <- result_output
    contrast_results[[collection_id]] <- result_output

    significant_result <- result_output[
      result_output$significant_fdr_0_05,
      ,
      drop = FALSE
    ]
    if (nrow(significant_result) > 0L) {
      for (row_index in seq_len(nrow(significant_result))) {
        entrez_ids <- strsplit(
          significant_result$leading_edge_entrez_ids[row_index],
          ";",
          fixed = TRUE
        )[[1L]]
        entrez_ids <- entrez_ids[entrez_ids != ""]
        if (length(entrez_ids) == 0L) {
          next
        }
        leading_edge_index <- leading_edge_index + 1L
        all_leading_edge[[leading_edge_index]] <- data.frame(
          contrast = significant_result$contrast[row_index],
          collection_id =
            significant_result$collection_id[row_index],
          collection_name =
            significant_result$collection_name[row_index],
          pathway = significant_result$pathway[row_index],
          NES =
            significant_result$normalized_enrichment_score[row_index],
          p_value = significant_result$p_value[row_index],
          FDR = significant_result$FDR[row_index],
          gene_entrez_id = entrez_ids,
          gene_symbol = unname(symbol_by_entrez[entrez_ids]),
          stringsAsFactors = FALSE
        )
      }
    }

    summary_index <- summary_index + 1L
    collection_summaries[[summary_index]] <- data.frame(
      contrast = contrast,
      collection_id = collection_id,
      collection_name = collection_name,
      source_gene_sets = length(resource$pathways),
      tested_gene_sets = nrow(result_output),
      minimum_tested_size = min(result_output$measured_gene_count),
      maximum_tested_size = max(result_output$measured_gene_count),
      significant_fdr_0_05 = nrow(significant_result),
      significant_positive_nes = sum(
        significant_result$normalized_enrichment_score > 0
      ),
      significant_negative_nes = sum(
        significant_result$normalized_enrichment_score < 0
      ),
      stringsAsFactors = FALSE
    )
  }

  contrast_output <- do.call(rbind, contrast_results)
  row.names(contrast_output) <- NULL
  write.csv(
    contrast_output,
    file.path(
      contrast_output_dir,
      paste0(contrast, "_gsea_results.csv")
    ),
    row.names = FALSE,
    na = ""
  )
}

combined_results <- do.call(rbind, all_results)
row.names(combined_results) <- NULL
significant_results <- combined_results[
  combined_results$significant_fdr_0_05,
  ,
  drop = FALSE
]
collection_summary <- do.call(rbind, collection_summaries)
row.names(collection_summary) <- NULL

if (length(all_leading_edge) == 0L) {
  leading_edge_genes <- make_empty_leading_edge()
} else {
  leading_edge_genes <- do.call(rbind, all_leading_edge)
  leading_edge_genes <- unique(leading_edge_genes)
  leading_edge_genes <- leading_edge_genes[
    order(
      leading_edge_genes$contrast,
      leading_edge_genes$collection_id,
      leading_edge_genes$FDR,
      leading_edge_genes$pathway,
      leading_edge_genes$gene_entrez_id
    ),
    ,
    drop = FALSE
  ]
  row.names(leading_edge_genes) <- NULL
}

write.csv(
  combined_results,
  combined_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  significant_results,
  significant_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  leading_edge_genes,
  leading_edge_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  collection_summary,
  collection_summary_file,
  row.names = FALSE,
  na = ""
)

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

if (nrow(significant_results) > 0L) {
  plot_split <- split(
    significant_results,
    interaction(
      significant_results$contrast,
      significant_results$collection_id,
      drop = TRUE
    )
  )
  plot_data <- do.call(
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
        head(data, 10L)
      }
    )
  )
  plot_data$pathway_label <- shorten_pathway_name(
    plot_data$pathway
  )
  plot_data$plot_id <- paste(
    plot_data$contrast,
    plot_data$collection_id,
    plot_data$pathway,
    sep = "|||"
  )
  plot_data <- plot_data[
    order(
      plot_data$collection_name,
      plot_data$contrast,
      plot_data$normalized_enrichment_score
    ),
    ,
    drop = FALSE
  ]
  plot_levels <- unique(plot_data$plot_id)
  plot_labels <- setNames(
    plot_data$pathway_label,
    plot_data$plot_id
  )
  plot_data$plot_id <- factor(
    plot_data$plot_id,
    levels = plot_levels
  )
  plot_data$minus_log10_fdr <- -log10(
    pmax(plot_data$FDR, .Machine$double.xmin)
  )
  plot_data$direction <- ifelse(
    plot_data$normalized_enrichment_score > 0,
    "Positive NES",
    "Negative NES"
  )

  summary_plot <- ggplot2::ggplot(
    plot_data,
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
    ggplot2::scale_y_discrete(labels = plot_labels) +
    ggplot2::scale_color_manual(
      values = c(
        "Negative NES" = "#2166AC",
        "Positive NES" = "#B2182B"
      )
    ) +
    ggplot2::labs(
      title = "FDR-significant mRNA GSEA pathways",
      subtitle = paste(
        "Up to 10 pathways with the smallest FDR per",
        "contrast and collection"
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
} else {
  summary_plot <- ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0,
      y = 0,
      label = "No pathways met FDR < 0.05."
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = "mRNA GSEA pathway summary") +
    ggplot2::theme_void()
}
ggplot2::ggsave(
  filename = summary_figure_file,
  plot = summary_plot,
  width = 16,
  height = 12,
  units = "in",
  dpi = 300,
  bg = "white"
)

summary_lines <- c(
  "mRNA GSEA run summary",
  "",
  paste0("Run date: ", Sys.Date()),
  paste0("mRNA DE model: treatment + sex + age"),
  paste0("Contrasts: ", paste(contrast_order, collapse = ", ")),
  paste0("Ranked genes per contrast: 20,082 unique Entrez genes"),
  paste0("Rank metric: limma moderated t-statistic"),
  paste(
    "Exact-tie rule: retain Script 11 rank order and apply only the",
    "smallest deterministic decrement needed to distinguish tied genes"
  ),
  "",
  paste0("GSEA method: fgsea::fgseaMultilevel"),
  paste0("fgsea version: ", as.character(packageVersion("fgsea"))),
  paste0("Minimum measured gene-set size: ", minimum_gene_set_size),
  paste0("Maximum measured gene-set size: ", maximum_gene_set_size),
  paste0("Significance threshold: FDR < ", fdr_threshold),
  paste0("Random seed base: ", random_seed),
  paste0("Parallel workers: 1"),
  "",
  paste0("MSigDB release: ", msigdb_version),
  paste0(
    "FDR scope: calculated separately for each collection within ",
    "each contrast"
  ),
  "",
  "Exact rank ties by contrast:"
)
tie_summary <- do.call(rbind, tie_summaries)
for (row_index in seq_len(nrow(tie_summary))) {
  row <- tie_summary[row_index, , drop = FALSE]
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ",
      row$contrast,
      ": ",
      row$exact_tie_groups,
      " exact tie groups containing ",
      row$genes_in_exact_ties,
      " genes; maximum numerical adjustment = ",
      format(
        row$maximum_absolute_tie_adjustment,
        scientific = TRUE,
        digits = 4
      )
    )
  )
}
summary_lines <- c(summary_lines, "", "Pathway collections:")
for (
  collection_index in seq_len(nrow(collection_definitions))
) {
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ",
      collection_definitions$collection_name[collection_index],
      " [",
      collection_definitions$collection_id[collection_index],
      "]: ",
      basename(collection_definitions$gmt_file[collection_index]),
      "; SHA-256 ",
      collection_definitions$sha256[collection_index]
    )
  )
}
summary_lines <- c(summary_lines, "", "Results by contrast and collection:")
for (row_index in seq_len(nrow(collection_summary))) {
  row <- collection_summary[row_index, , drop = FALSE]
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ",
      row$contrast,
      " / ",
      row$collection_name,
      ": ",
      row$tested_gene_sets,
      " tested; ",
      row$significant_fdr_0_05,
      " FDR-significant (",
      row$significant_positive_nes,
      " positive NES, ",
      row$significant_negative_nes,
      " negative NES)"
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  paste0("Combined tested pathway rows: ", nrow(combined_results)),
  paste0(
    "FDR-significant pathway rows: ",
    nrow(significant_results)
  ),
  paste0(
    "Leading-edge gene rows from FDR-significant pathways: ",
    nrow(leading_edge_genes)
  ),
  "",
  paste0("Wrote contrast results under: ", contrast_output_dir),
  paste0("Wrote combined results: ", combined_output_file),
  paste0("Wrote significant results: ", significant_output_file),
  paste0("Wrote leading-edge genes: ", leading_edge_output_file),
  paste0("Wrote collection summary: ", collection_summary_file),
  paste0("Wrote summary figure: ", summary_figure_file),
  "",
  paste(
    "Interpretation boundary: GSEA identifies coordinated mRNA",
    "pathway-level expression patterns. It does not establish direct",
    "miRNA regulation, causality, or mediation."
  ),
  paste0("R version: ", R.version.string)
)
writeLines(summary_lines, summary_output_file)

cat(
  "GSEA tested ",
  nrow(combined_results),
  " contrast--collection pathway rows; ",
  nrow(significant_results),
  " met FDR < ",
  fdr_threshold,
  ".\n",
  sep = ""
)
cat("Wrote: ", summary_output_file, "\n", sep = "")
