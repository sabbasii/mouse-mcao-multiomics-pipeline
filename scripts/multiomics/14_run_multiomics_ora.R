#!/usr/bin/env Rscript

# Run formal over-representation analysis (ORA) of the 53 selected genes.
#
# This script performs ORA only. It does not run GSEA and does not use a
# genome-wide background. The tested list contains 53 unique mRNA genes from
# the 56 DE-supported miRNA--mRNA pairs in MCAO3hr_vs_MCAO1hr. The statistical
# background contains 1,106 unique measured genes targeted by at least one of
# the same six selected miRNAs. This background reflects the genes that could
# have entered the selected list under the target-linking design.
#
# Mouse MSigDB 2025.1 Hallmark and canonical-pathway collections are tested
# separately. A pathway is eligible when 5 through 500 of its genes occur in
# the 1,106-gene background. One-sided Fisher exact tests ask whether selected
# genes are over-represented in each pathway. Benjamini-Hochberg FDR correction
# is applied separately to each collection.
#
# Inputs:
#   results/multiomics/enrichment_inputs/ora/
#     formal_ora_selected_53_target_linked_genes.csv
#   results/multiomics/enrichment_inputs/backgrounds/
#     03_candidate_mirna_target_gene_background.csv
#   resources/pathway_gene_sets/msigdb_2025.1.Mm/
#     mh.all.v2025.1.Mm.entrez.gmt
#     m2.cp.v2025.1.Mm.entrez.gmt
#
# Outputs:
#   results/multiomics/enrichment_results/ora/
#     collection_results/mouse_hallmark_ora_results.csv
#     collection_results/mouse_canonical_pathways_ora_results.csv
#     combined_ora_results.csv
#     significant_ora_results_fdr_0.05.csv
#     nominal_p_lt_0.05_ora_results.csv
#     pathway_selected_gene_overlaps.csv
#     pathway_collection_summary.csv
#     figures/ora_top_pathways.png
#     ora_run_summary.txt
#
# Usage:
#   Rscript scripts/multiomics/14_run_multiomics_ora.R

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

selected_gene_file <- file.path(
  "results", "multiomics", "enrichment_inputs", "ora",
  "formal_ora_selected_53_target_linked_genes.csv"
)
background_gene_file <- file.path(
  "results", "multiomics", "enrichment_inputs", "backgrounds",
  "03_candidate_mirna_target_gene_background.csv"
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
  output_stem = c(
    "mouse_hallmark",
    "mouse_canonical_pathways"
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
  "results", "multiomics", "enrichment_results", "ora"
)
collection_output_dir <- file.path(output_dir, "collection_results")
figure_output_dir <- file.path(output_dir, "figures")
combined_output_file <- file.path(
  output_dir, "combined_ora_results.csv"
)
significant_output_file <- file.path(
  output_dir, "significant_ora_results_fdr_0.05.csv"
)
nominal_output_file <- file.path(
  output_dir, "nominal_p_lt_0.05_ora_results.csv"
)
overlap_output_file <- file.path(
  output_dir, "pathway_selected_gene_overlaps.csv"
)
collection_summary_file <- file.path(
  output_dir, "pathway_collection_summary.csv"
)
figure_output_file <- file.path(
  figure_output_dir, "ora_top_pathways.png"
)
run_summary_output_file <- file.path(
  output_dir, "ora_run_summary.txt"
)

minimum_background_pathway_size <- 5L
maximum_background_pathway_size <- 500L
fdr_threshold <- 0.05
nominal_p_threshold <- 0.05

required_inputs <- c(
  selected_gene_file,
  background_gene_file,
  collection_definitions$gmt_file
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

dir.create(
  collection_output_dir,
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
  if (any(lengths(fields) < 3L)) {
    stop("Malformed GMT line(s) in: ", path, call. = FALSE)
  }
  pathway_names <- vapply(fields, `[[`, character(1), 1L)
  descriptions <- vapply(fields, `[[`, character(1), 2L)
  if (
    any(pathway_names == "") ||
      anyDuplicated(pathway_names)
  ) {
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
  metadata <- data.frame(
    pathway = pathway_names,
    pathway_description = descriptions,
    source_gene_count = lengths(pathways),
    stringsAsFactors = FALSE
  )
  list(pathways = pathways, metadata = metadata)
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

selected_genes <- read.csv(
  selected_gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
background_genes <- read.csv(
  background_gene_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_selected_columns <- c(
  "gene_entrez_id", "gene_symbol", "selected_mirna_ids",
  "mrna_logFC", "mrna_P.Value", "mrna_adj.P.Val"
)
required_background_columns <- c(
  "gene_entrez_id", "gene_symbol", "measured_mirna_ids"
)
if (
  length(setdiff(
    required_selected_columns,
    names(selected_genes)
  )) > 0L ||
    length(setdiff(
      required_background_columns,
      names(background_genes)
    )) > 0L
) {
  stop("ORA gene input columns failed validation.", call. = FALSE)
}
selected_genes$gene_entrez_id <- trimws(
  as.character(selected_genes$gene_entrez_id)
)
background_genes$gene_entrez_id <- trimws(
  as.character(background_genes$gene_entrez_id)
)
split_mirna_ids <- function(x) {
  unique(unlist(strsplit(x, " || ", fixed = TRUE)))
}
selected_mirnas <- sort(split_mirna_ids(
  selected_genes$selected_mirna_ids
))
background_mirnas <- sort(split_mirna_ids(
  background_genes$measured_mirna_ids
))
if (
  nrow(selected_genes) != 53L ||
    nrow(background_genes) != 1106L ||
    anyNA(selected_genes$gene_entrez_id) ||
    anyNA(background_genes$gene_entrez_id) ||
    anyDuplicated(selected_genes$gene_entrez_id) ||
    anyDuplicated(background_genes$gene_entrez_id) ||
    length(selected_mirnas) != 6L ||
    !identical(selected_mirnas, background_mirnas) ||
    !all(
      selected_genes$gene_entrez_id %in%
        background_genes$gene_entrez_id
    )
) {
  stop(
    "The 53-gene ORA list or 1,106-gene background failed validation.",
    call. = FALSE
  )
}

selected_ids <- selected_genes$gene_entrez_id
background_ids <- background_genes$gene_entrez_id
selected_size <- length(selected_ids)
background_size <- length(background_ids)
background_nonselected_size <- background_size - selected_size
symbol_by_entrez <- setNames(
  background_genes$gene_symbol,
  background_genes$gene_entrez_id
)
mirnas_by_entrez <- setNames(
  background_genes$measured_mirna_ids,
  background_genes$gene_entrez_id
)

collection_results <- list()
overlap_tables <- list()
collection_summaries <- list()

for (collection_index in seq_len(nrow(collection_definitions))) {
  collection_id <- collection_definitions$collection_id[
    collection_index
  ]
  collection_name <- collection_definitions$collection_name[
    collection_index
  ]
  resource <- read_gmt(
    collection_definitions$gmt_file[collection_index]
  )
  background_pathway_genes <- lapply(
    resource$pathways,
    intersect,
    background_ids
  )
  background_pathway_sizes <- lengths(background_pathway_genes)
  eligible <- background_pathway_sizes >=
    minimum_background_pathway_size &
    background_pathway_sizes <= maximum_background_pathway_size
  eligible_pathways <- names(background_pathway_genes)[eligible]
  if (length(eligible_pathways) == 0L) {
    stop(
      "No eligible pathways for collection: ",
      collection_name,
      call. = FALSE
    )
  }

  result_rows <- lapply(
    eligible_pathways,
    function(pathway) {
      pathway_background_ids <- background_pathway_genes[[pathway]]
      pathway_selected_ids <- intersect(
        pathway_background_ids,
        selected_ids
      )
      selected_in_pathway <- length(pathway_selected_ids)
      selected_not_in_pathway <-
        selected_size - selected_in_pathway
      nonselected_in_pathway <-
        length(pathway_background_ids) - selected_in_pathway
      nonselected_not_in_pathway <-
        background_nonselected_size - nonselected_in_pathway
      contingency_table <- matrix(
        c(
          selected_in_pathway,
          selected_not_in_pathway,
          nonselected_in_pathway,
          nonselected_not_in_pathway
        ),
        nrow = 2L,
        byrow = TRUE,
        dimnames = list(
          selection = c("selected", "not_selected"),
          pathway = c("in_pathway", "not_in_pathway")
        )
      )
      fisher_result <- fisher.test(
        contingency_table,
        alternative = "greater"
      )
      background_pathway_size <- length(pathway_background_ids)
      selected_gene_ratio <-
        selected_in_pathway / selected_size
      background_gene_ratio <-
        background_pathway_size / background_size
      fold_enrichment <- if (background_gene_ratio == 0) {
        NA_real_
      } else {
        selected_gene_ratio / background_gene_ratio
      }
      metadata_index <- match(pathway, resource$metadata$pathway)
      data.frame(
        collection_id = collection_id,
        collection_name = collection_name,
        msigdb_version = msigdb_version,
        pathway = pathway,
        pathway_description =
          resource$metadata$pathway_description[metadata_index],
        source_gene_count =
          resource$metadata$source_gene_count[metadata_index],
        background_gene_count = background_pathway_size,
        selected_gene_count = selected_in_pathway,
        selected_gene_ratio = selected_gene_ratio,
        background_gene_ratio = background_gene_ratio,
        fold_enrichment = fold_enrichment,
        odds_ratio = unname(fisher_result$estimate),
        p_value = fisher_result$p.value,
        selected_gene_entrez_ids = paste(
          pathway_selected_ids,
          collapse = ";"
        ),
        selected_gene_symbols = paste(
          unname(symbol_by_entrez[pathway_selected_ids]),
          collapse = ";"
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  result <- do.call(rbind, result_rows)
  result$FDR <- p.adjust(result$p_value, method = "BH")
  result$nominal_p_lt_0_05 <-
    result$p_value < nominal_p_threshold
  result$significant_fdr_0_05 <-
    result$FDR < fdr_threshold
  result <- result[
    order(
      result$FDR,
      result$p_value,
      -result$selected_gene_count,
      -result$fold_enrichment,
      result$pathway
    ),
    ,
    drop = FALSE
  ]
  row.names(result) <- NULL
  collection_results[[collection_id]] <- result

  overlap_result <- result[
    result$selected_gene_count > 0L,
    ,
    drop = FALSE
  ]
  overlap_rows <- list()
  overlap_index <- 0L
  if (nrow(overlap_result) > 0L) {
    for (row_index in seq_len(nrow(overlap_result))) {
      pathway_selected_ids <- strsplit(
        overlap_result$selected_gene_entrez_ids[row_index],
        ";",
        fixed = TRUE
      )[[1L]]
      pathway_selected_ids <- pathway_selected_ids[
        pathway_selected_ids != ""
      ]
      for (gene_id in pathway_selected_ids) {
        selected_index <- match(
          gene_id,
          selected_genes$gene_entrez_id
        )
        overlap_index <- overlap_index + 1L
        overlap_rows[[overlap_index]] <- data.frame(
          collection_id = collection_id,
          collection_name = collection_name,
          pathway = overlap_result$pathway[row_index],
          pathway_p_value =
            overlap_result$p_value[row_index],
          pathway_FDR = overlap_result$FDR[row_index],
          pathway_fold_enrichment =
            overlap_result$fold_enrichment[row_index],
          gene_entrez_id = gene_id,
          gene_symbol = unname(symbol_by_entrez[gene_id]),
          selected_mirna_ids =
            selected_genes$selected_mirna_ids[selected_index],
          mrna_logFC =
            selected_genes$mrna_logFC[selected_index],
          mrna_P.Value =
            selected_genes$mrna_P.Value[selected_index],
          mrna_adj.P.Val =
            selected_genes$mrna_adj.P.Val[selected_index],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(overlap_rows) == 0L) {
    overlap_tables[[collection_id]] <- data.frame(
      collection_id = character(0),
      collection_name = character(0),
      pathway = character(0),
      pathway_p_value = numeric(0),
      pathway_FDR = numeric(0),
      pathway_fold_enrichment = numeric(0),
      gene_entrez_id = character(0),
      gene_symbol = character(0),
      selected_mirna_ids = character(0),
      mrna_logFC = numeric(0),
      mrna_P.Value = numeric(0),
      mrna_adj.P.Val = numeric(0),
      stringsAsFactors = FALSE
    )
  } else {
    overlap_tables[[collection_id]] <- do.call(
      rbind,
      overlap_rows
    )
  }

  collection_summaries[[collection_id]] <- data.frame(
    collection_id = collection_id,
    collection_name = collection_name,
    source_pathway_count = length(resource$pathways),
    tested_pathway_count = nrow(result),
    pathways_with_selected_gene_overlap = sum(
      result$selected_gene_count > 0L
    ),
    nominal_p_lt_0_05 = sum(result$nominal_p_lt_0_05),
    significant_fdr_0_05 = sum(result$significant_fdr_0_05),
    stringsAsFactors = FALSE
  )

  write.csv(
    result,
    file.path(
      collection_output_dir,
      paste0(
        collection_definitions$output_stem[collection_index],
        "_ora_results.csv"
      )
    ),
    row.names = FALSE,
    na = ""
  )
}

combined_results <- do.call(rbind, collection_results)
row.names(combined_results) <- NULL
significant_results <- combined_results[
  combined_results$significant_fdr_0_05,
  ,
  drop = FALSE
]
nominal_results <- combined_results[
  combined_results$nominal_p_lt_0_05,
  ,
  drop = FALSE
]
overlap_results <- do.call(rbind, overlap_tables)
row.names(overlap_results) <- NULL
collection_summary <- do.call(rbind, collection_summaries)
row.names(collection_summary) <- NULL

if (
  nrow(combined_results) !=
    sum(collection_summary$tested_pathway_count) ||
    anyNA(combined_results$p_value) ||
    anyNA(combined_results$FDR) ||
    any(combined_results$p_value < 0 |
      combined_results$p_value > 1) ||
    any(combined_results$FDR < 0 |
      combined_results$FDR > 1) ||
    any(
      combined_results$selected_gene_count >
        combined_results$background_gene_count
    ) ||
    anyDuplicated(
      paste(
        combined_results$collection_id,
        combined_results$pathway
      )
    )
) {
  stop("ORA result validation failed.", call. = FALSE)
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
  nominal_results,
  nominal_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  overlap_results,
  overlap_output_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  collection_summary,
  collection_summary_file,
  row.names = FALSE,
  na = ""
)

plot_split <- split(
  combined_results,
  combined_results$collection_id
)
plot_data <- do.call(
  rbind,
  lapply(
    plot_split,
    function(data) {
      data <- data[
        order(
          data$p_value,
          data$FDR,
          -data$selected_gene_count,
          -data$fold_enrichment
        ),
        ,
        drop = FALSE
      ]
      head(data, 10L)
    }
  )
)
plot_data$pathway_label <- shorten_pathway_name(plot_data$pathway)
plot_data$plot_id <- paste(
  plot_data$collection_id,
  plot_data$pathway,
  sep = "|||"
)
plot_data <- plot_data[
  order(
    plot_data$collection_name,
    plot_data$fold_enrichment
  ),
  ,
  drop = FALSE
]
plot_data$plot_id <- factor(
  plot_data$plot_id,
  levels = unique(plot_data$plot_id)
)
plot_labels <- setNames(
  plot_data$pathway_label,
  plot_data$plot_id
)
plot_data$nominal_status <- ifelse(
  plot_data$nominal_p_lt_0_05,
  "Raw p < 0.05",
  "Raw p >= 0.05"
)

plot_subtitle <- if (nrow(significant_results) == 0L) {
  paste(
    "Top 10 pathways by raw p-value per collection;",
    "none met FDR < 0.05"
  )
} else {
  paste(
    "Top 10 pathways by raw p-value per collection;",
    "FDR-significant pathways are present"
  )
}
ora_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = fold_enrichment,
    y = plot_id,
    color = nominal_status,
    size = selected_gene_count
  )
) +
  ggplot2::geom_vline(
    xintercept = 1,
    color = "grey60",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  ggplot2::geom_point(alpha = 0.85) +
  ggplot2::facet_wrap(
    ggplot2::vars(collection_name),
    ncol = 2,
    scales = "free_y"
  ) +
  ggplot2::scale_y_discrete(labels = plot_labels) +
  ggplot2::scale_color_manual(
    values = c(
      "Raw p < 0.05" = "#B2182B",
      "Raw p >= 0.05" = "#636363"
    )
  ) +
  ggplot2::labs(
    title = "Formal ORA of 53 selected target-linked genes",
    subtitle = plot_subtitle,
    x = "Fold enrichment relative to the 1,106-gene background",
    y = NULL,
    color = NULL,
    size = "Selected genes\nin pathway"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )
ggplot2::ggsave(
  filename = figure_output_file,
  plot = ora_plot,
  width = 14,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

summary_lines <- c(
  "Formal multi-omics ORA run summary",
  "",
  paste0("Run date: ", Sys.Date()),
  paste0("Selected contrast: MCAO3hr_vs_MCAO1hr"),
  paste0(
    "Selection rule: miRNA FDR < 0.10 and mRNA FDR < 0.10, ",
    "plus measured miRNA--target evidence"
  ),
  paste0("Selected genes: ", selected_size),
  paste0("Background genes: ", background_size),
  paste0(
    "Background definition: measured targets of the same six ",
    "selected miRNAs"
  ),
  "",
  paste0("Test: one-sided Fisher exact over-representation test"),
  paste0(
    "Eligible background pathway size: ",
    minimum_background_pathway_size,
    " through ",
    maximum_background_pathway_size
  ),
  paste0(
    "FDR scope: Benjamini-Hochberg correction separately within ",
    "each MSigDB collection"
  ),
  paste0("Formal significance threshold: FDR < ", fdr_threshold),
  paste0(
    "Raw p < ",
    nominal_p_threshold,
    " is reported descriptively and is not formal significance"
  ),
  "",
  paste0("MSigDB release: ", msigdb_version)
)
for (collection_index in seq_len(nrow(collection_definitions))) {
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ",
      collection_definitions$collection_name[collection_index],
      ": ",
      basename(collection_definitions$gmt_file[collection_index]),
      "; SHA-256 ",
      collection_definitions$sha256[collection_index]
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "Results by collection:"
)
for (row_index in seq_len(nrow(collection_summary))) {
  row <- collection_summary[row_index, , drop = FALSE]
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ",
      row$collection_name,
      ": ",
      row$tested_pathway_count,
      " tested; ",
      row$nominal_p_lt_0_05,
      " raw p < 0.05; ",
      row$significant_fdr_0_05,
      " FDR < 0.05"
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  paste0("Combined tested pathway rows: ", nrow(combined_results)),
  paste0("Raw p < 0.05 pathway rows: ", nrow(nominal_results)),
  paste0(
    "FDR < 0.05 pathway rows: ",
    nrow(significant_results)
  ),
  paste0("Selected pathway--gene overlap rows: ", nrow(overlap_results)),
  "",
  paste0("Wrote collection results under: ", collection_output_dir),
  paste0("Wrote combined results: ", combined_output_file),
  paste0("Wrote FDR-significant results: ", significant_output_file),
  paste0("Wrote nominal results: ", nominal_output_file),
  paste0("Wrote pathway--gene overlaps: ", overlap_output_file),
  paste0("Wrote collection summary: ", collection_summary_file),
  paste0("Wrote summary figure: ", figure_output_file),
  "",
  paste(
    "Interpretation boundary: ORA tests whether the 53 selected genes",
    "occur in pathways more often than expected among the 1,106 genes",
    "that could have been selected. A raw p-value below 0.05 without",
    "FDR support is exploratory and not a statistically significant",
    "pathway result."
  ),
  paste(
    "ORA does not establish direct miRNA regulation, causality, or",
    "mediation."
  ),
  paste0("R version: ", R.version.string)
)
writeLines(summary_lines, run_summary_output_file)

cat(
  "ORA tested ",
  nrow(combined_results),
  " pathways; ",
  nrow(significant_results),
  " met FDR < ",
  fdr_threshold,
  ".\n",
  sep = ""
)
cat("Wrote: ", run_summary_output_file, "\n", sep = "")
