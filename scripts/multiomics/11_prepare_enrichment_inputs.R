#!/usr/bin/env Rscript

# Prepare gene-level inputs for subsequent enrichment analyses.
#
# This script prepares, but does not run:
#   - conventional mRNA GSEA using all measured genes ranked by moderated t;
#   - formal over-representation analysis (ORA) of the 53 selected genes;
#   - miRNA-target gene-set enrichment for the six selected miRNAs.
#
# The 12 inverse-DE genes and 9 directionally compatible but unconfirmed
# genes are retained as descriptive review lists only. They are too small for
# reliable standalone ORA and must not be treated as formal enrichment inputs.
#
# The mRNA DE tables contain transcript-cluster-level results. When more than
# one uniquely mapped transcript cluster represents the same Entrez gene, one
# feature is selected using the highest average expression. This rule is
# applied once and does not use the contrast-specific t-statistic or p-value.
#
# Inputs:
#   results/mrna/differential_expression/treatment_sex_age_limma/
#     transcript_cluster_<contrast>_mrna.csv
#   results/multiomics/miRNA_target_evidence/analysis_ready/
#     mirna_mrna_unique_pair_summary.csv
#   results/multiomics/integrated_pair_evidence/
#     mirna_mrna_integrated_pair_evidence_summary.csv
#     mirna_mrna_directionally_compatible_unconfirmed_pairs.csv
#
# Outputs:
#   results/multiomics/enrichment_inputs/
#     gsea/
#     ora/
#     descriptive_review_lists/
#     backgrounds/
#     mirna_target_gene_sets/
#     enrichment_input_summary.txt
#
# Usage:
#   Rscript scripts/multiomics/11_prepare_enrichment_inputs.R

mrna_de_dir <- file.path(
  "results", "mrna", "differential_expression",
  "treatment_sex_age_limma"
)
unique_pair_file <- file.path(
  "results", "multiomics", "miRNA_target_evidence", "analysis_ready",
  "mirna_mrna_unique_pair_summary.csv"
)
integrated_pair_file <- file.path(
  "results", "multiomics", "integrated_pair_evidence",
  "mirna_mrna_integrated_pair_evidence_summary.csv"
)
compatible_pair_file <- file.path(
  "results", "multiomics", "integrated_pair_evidence",
  "mirna_mrna_directionally_compatible_unconfirmed_pairs.csv"
)

contrast_order <- c(
  "MCAO1hr_vs_Sham",
  "MCAO3hr_vs_Sham",
  "MCAO3hr_vs_MCAO1hr"
)
mrna_de_files <- setNames(
  file.path(
    mrna_de_dir,
    paste0(
      "transcript_cluster_",
      contrast_order,
      "_mrna.csv"
    )
  ),
  contrast_order
)

output_dir <- file.path(
  "results", "multiomics", "enrichment_inputs"
)
gsea_dir <- file.path(output_dir, "gsea")
ora_dir <- file.path(output_dir, "ora")
review_list_dir <- file.path(
  output_dir, "descriptive_review_lists"
)
background_dir <- file.path(output_dir, "backgrounds")
target_set_dir <- file.path(
  output_dir, "mirna_target_gene_sets"
)
summary_output_file <- file.path(
  output_dir, "enrichment_input_summary.txt"
)

all_measured_background_file <- file.path(
  background_dir, "01_all_measured_mrna_genes.csv"
)
all_target_linked_background_file <- file.path(
  background_dir, "02_all_measured_mirna_target_genes.csv"
)
candidate_target_background_file <- file.path(
  background_dir,
  "03_candidate_mirna_target_gene_background.csv"
)

selected_gene_file <- file.path(
  ora_dir, "formal_ora_selected_53_target_linked_genes.csv"
)
inverse_gene_file <- file.path(
  review_list_dir,
  "inverse_de_12_unique_genes_descriptive.csv"
)
compatible_gene_file <- file.path(
  review_list_dir,
  paste0(
    "directionally_compatible_unconfirmed_",
    "9_unique_genes_descriptive.csv"
  )
)
review_note_file <- file.path(
  review_list_dir, "README.txt"
)

target_set_long_file <- file.path(
  target_set_dir,
  "candidate_mirna_measured_target_gene_sets_long.csv"
)
target_set_gmt_file <- file.path(
  target_set_dir,
  "candidate_mirna_measured_target_gene_sets_entrez.gmt"
)
target_set_size_file <- file.path(
  target_set_dir,
  "candidate_mirna_measured_target_gene_set_sizes.csv"
)

required_inputs <- c(
  unname(mrna_de_files),
  unique_pair_file,
  integrated_pair_file,
  compatible_pair_file
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

clean_id <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "null")] <- NA_character_
  x
}

collapse_values <- function(x, separator = " || ") {
  x <- trimws(as.character(x))
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x) == 0L) {
    return(NA_character_)
  }
  paste(x, collapse = separator)
}

mrna_de <- lapply(mrna_de_files, read_table)
required_de_columns <- c(
  "transcript_cluster_id", "logFC", "AveExpr", "t",
  "P.Value", "adj.P.Val", "B", "ENTREZID",
  "SYMBOL", "GENENAME", "GENETYPE", "ENSEMBL"
)
for (contrast in contrast_order) {
  data <- mrna_de[[contrast]]
  missing_columns <- setdiff(required_de_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      contrast,
      " mRNA DE table lacks required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (
    nrow(data) != 20222L ||
      anyNA(data[
        ,
        setdiff(required_de_columns, "ENSEMBL"),
        drop = FALSE
      ]) ||
      anyDuplicated(data$transcript_cluster_id) ||
      any(!is.finite(data$logFC)) ||
      any(!is.finite(data$AveExpr)) ||
      any(!is.finite(data$t)) ||
      any(!is.finite(data$P.Value)) ||
      any(!is.finite(data$adj.P.Val))
  ) {
    stop(
      contrast,
      " mRNA DE table failed completeness or uniqueness checks.",
      call. = FALSE
    )
  }
  data$ENTREZID <- clean_id(data$ENTREZID)
  data$SYMBOL <- clean_id(data$SYMBOL)
  data$ENSEMBL <- clean_id(data$ENSEMBL)
  if (
    anyNA(data$ENTREZID) ||
      anyNA(data$SYMBOL)
  ) {
    stop(
      contrast,
      " mRNA DE table has missing Entrez IDs or gene symbols.",
      call. = FALSE
    )
  }
  mrna_de[[contrast]] <- data
}

reference_de <- mrna_de[[contrast_order[1L]]]
reference_feature_order <- order(reference_de$transcript_cluster_id)
reference_de <- reference_de[
  reference_feature_order,
  ,
  drop = FALSE
]

constant_de_columns <- c(
  "transcript_cluster_id", "AveExpr", "ENTREZID",
  "SYMBOL", "GENENAME", "GENETYPE", "ENSEMBL"
)
for (contrast in contrast_order[-1L]) {
  data <- mrna_de[[contrast]]
  data <- data[
    match(reference_de$transcript_cluster_id, data$transcript_cluster_id),
    ,
    drop = FALSE
  ]
  for (column in constant_de_columns) {
    reference_value <- as.character(reference_de[[column]])
    comparison_value <- as.character(data[[column]])
    if (!identical(reference_value, comparison_value)) {
      stop(
        column,
        " differs across mRNA contrast tables.",
        call. = FALSE
      )
    }
  }
}

features_per_entrez <- table(reference_de$ENTREZID)
feature_selection_order <- order(
  reference_de$ENTREZID,
  -reference_de$AveExpr,
  reference_de$transcript_cluster_id
)
feature_selection <- reference_de[
  feature_selection_order,
  c(
    "transcript_cluster_id", "AveExpr", "ENTREZID",
    "SYMBOL", "GENENAME", "GENETYPE", "ENSEMBL"
  ),
  drop = FALSE
]
feature_selection$features_mapped_to_entrez <- as.integer(
  features_per_entrez[feature_selection$ENTREZID]
)
feature_selection$selected_for_gene_level_enrichment <-
  !duplicated(feature_selection$ENTREZID)
feature_selection$selection_rule <- ifelse(
  feature_selection$selected_for_gene_level_enrichment,
  "highest_AveExpr_then_transcript_cluster_id",
  "not_selected_lower_AveExpr_or_later_tie"
)

selected_features <- feature_selection[
  feature_selection$selected_for_gene_level_enrichment,
  ,
  drop = FALSE
]
selected_transcript_clusters <- selected_features$transcript_cluster_id
if (
  nrow(selected_features) != 20082L ||
    anyDuplicated(selected_features$ENTREZID) ||
    length(selected_transcript_clusters) != 20082L
) {
  stop(
    "Gene-level mRNA feature selection did not produce 20,082 genes.",
    call. = FALSE
  )
}

all_measured_background <- selected_features[
  ,
  c(
    "ENTREZID", "SYMBOL", "GENENAME", "GENETYPE",
    "ENSEMBL", "transcript_cluster_id", "AveExpr",
    "features_mapped_to_entrez", "selection_rule"
  ),
  drop = FALSE
]
names(all_measured_background)[
  names(all_measured_background) == "ENTREZID"
] <- "gene_entrez_id"
names(all_measured_background)[
  names(all_measured_background) == "SYMBOL"
] <- "gene_symbol"
names(all_measured_background)[
  names(all_measured_background) == "GENENAME"
] <- "gene_name"
names(all_measured_background)[
  names(all_measured_background) == "GENETYPE"
] <- "gene_type"
names(all_measured_background)[
  names(all_measured_background) == "ENSEMBL"
] <- "gene_ensembl_id"
names(all_measured_background)[
  names(all_measured_background) == "AveExpr"
] <- "average_expression"
all_measured_background <- all_measured_background[
  order(all_measured_background$gene_entrez_id),
  ,
  drop = FALSE
]
row.names(all_measured_background) <- NULL

gsea_output_files <- character(0)
gsea_tie_counts <- integer(0)
for (contrast in contrast_order) {
  data <- mrna_de[[contrast]]
  selected_index <- match(
    selected_transcript_clusters,
    data$transcript_cluster_id
  )
  if (anyNA(selected_index)) {
    stop(
      "Selected transcript clusters are missing from ", contrast, ".",
      call. = FALSE
    )
  }

  ranked <- data[
    selected_index,
    required_de_columns,
    drop = FALSE
  ]
  ranked$rank_metric <- ranked$t
  ranked <- ranked[
    order(-ranked$rank_metric, ranked$ENTREZID),
    ,
    drop = FALSE
  ]
  ranked$rank_position <- seq_len(nrow(ranked))
  if (
    nrow(ranked) != 20082L ||
      anyDuplicated(ranked$ENTREZID) ||
      any(!is.finite(ranked$rank_metric))
  ) {
    stop(
      contrast,
      " GSEA ranking is incomplete or contains duplicated genes.",
      call. = FALSE
    )
  }

  csv_file <- file.path(
    gsea_dir,
    paste0("gsea_ranked_entrez_", contrast, ".csv")
  )
  rnk_file <- file.path(
    gsea_dir,
    paste0("gsea_ranked_entrez_", contrast, ".rnk")
  )
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(ranked, csv_file, row.names = FALSE, na = "")
  write.table(
    ranked[, c("ENTREZID", "rank_metric"), drop = FALSE],
    rnk_file,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    na = ""
  )
  gsea_output_files <- c(gsea_output_files, csv_file, rnk_file)
  gsea_tie_counts[contrast] <- sum(duplicated(ranked$rank_metric))
}

unique_pairs <- read_table(unique_pair_file)
integrated_pairs <- read_table(integrated_pair_file)
compatible_pairs <- read_table(compatible_pair_file)

required_unique_pair_columns <- c(
  "mirna_id", "gene_entrez_id", "gene_symbol",
  "gene_name", "gene_type", "gene_ensembl_id",
  "mirtarbase_supported", "targetscanmouse_supported"
)
missing_unique_pair_columns <- setdiff(
  required_unique_pair_columns,
  names(unique_pairs)
)
if (length(missing_unique_pair_columns) > 0L) {
  stop(
    "Unique-pair table lacks required column(s): ",
    paste(missing_unique_pair_columns, collapse = ", "),
    call. = FALSE
  )
}

required_integrated_columns <- c(
  "mirna_id", "gene_entrez_id", "gene_symbol",
  "gene_name", "gene_type", "gene_ensembl_id",
  "mrna_logFC", "mrna_P.Value", "mrna_adj.P.Val",
  "inverse_logFC_direction",
  "directionally_compatible_unconfirmed",
  "target_evidence_category"
)
missing_integrated_columns <- setdiff(
  required_integrated_columns,
  names(integrated_pairs)
)
if (length(missing_integrated_columns) > 0L) {
  stop(
    "Integrated-pair table lacks required column(s): ",
    paste(missing_integrated_columns, collapse = ", "),
    call. = FALSE
  )
}

if (
  nrow(integrated_pairs) != 56L ||
    nrow(compatible_pairs) != 10L
) {
  stop(
    "Integrated pair inputs do not have the expected 56 and 10 rows.",
    call. = FALSE
  )
}

for (data_name in c(
  "unique_pairs", "integrated_pairs", "compatible_pairs"
)) {
  data <- get(data_name)
  data$gene_entrez_id <- clean_id(data$gene_entrez_id)
  data$gene_symbol <- clean_id(data$gene_symbol)
  data$gene_ensembl_id <- clean_id(data$gene_ensembl_id)
  if (
    anyNA(data$gene_entrez_id) ||
      anyNA(data$gene_symbol)
  ) {
    stop(
      data_name,
      " contains missing Entrez IDs or gene symbols.",
      call. = FALSE
    )
  }
  assign(data_name, data)
}

build_target_gene_background <- function(data, candidate_mirnas = NULL) {
  if (!is.null(candidate_mirnas)) {
    data <- data[data$mirna_id %in% candidate_mirnas, , drop = FALSE]
  }
  split_index <- split(
    seq_len(nrow(data)),
    data$gene_entrez_id
  )
  output <- lapply(split_index, function(index) {
    group <- data[index, , drop = FALSE]
    gene_symbols <- unique(group$gene_symbol)
    gene_names <- unique(group$gene_name)
    gene_types <- unique(group$gene_type)
    gene_ensembl_ids <- unique(
      group$gene_ensembl_id[!is.na(group$gene_ensembl_id)]
    )
    if (
      length(gene_symbols) != 1L ||
        length(gene_names) != 1L ||
        length(gene_types) != 1L ||
        length(gene_ensembl_ids) > 1L
    ) {
      stop(
        "Gene annotation varies within a target-linked Entrez gene.",
        call. = FALSE
      )
    }
    data.frame(
      gene_entrez_id = group$gene_entrez_id[1L],
      gene_symbol = gene_symbols,
      gene_name = gene_names,
      gene_type = gene_types,
      gene_ensembl_id = if (length(gene_ensembl_ids) == 0L) {
        NA_character_
      } else {
        gene_ensembl_ids
      },
      measured_mirna_count = length(unique(group$mirna_id)),
      measured_mirna_ids = collapse_values(group$mirna_id),
      mirtarbase_supported = any(group$mirtarbase_supported),
      targetscanmouse_supported = any(
        group$targetscanmouse_supported
      ),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, output)
  output$target_evidence_category <- ifelse(
    output$mirtarbase_supported &
      output$targetscanmouse_supported,
    "miRTarBase_and_TargetScanMouse",
    ifelse(
      output$mirtarbase_supported,
      "miRTarBase_only",
      "TargetScanMouse_only"
    )
  )
  output <- output[
    order(output$gene_entrez_id),
    ,
    drop = FALSE
  ]
  row.names(output) <- NULL
  output
}

all_target_linked_background <- build_target_gene_background(
  unique_pairs
)
candidate_mirnas <- sort(unique(integrated_pairs$mirna_id))
candidate_target_background <- build_target_gene_background(
  unique_pairs,
  candidate_mirnas = candidate_mirnas
)

if (
  nrow(all_target_linked_background) != 12742L ||
    nrow(candidate_target_background) != 1106L ||
    length(candidate_mirnas) != 6L ||
    !all(
      all_target_linked_background$gene_entrez_id %in%
        all_measured_background$gene_entrez_id
    ) ||
    !all(
      candidate_target_background$gene_entrez_id %in%
        all_target_linked_background$gene_entrez_id
    )
) {
  stop(
    "Measured target-linked gene backgrounds failed validation.",
    call. = FALSE
  )
}

collapse_selected_genes <- function(data) {
  split_index <- split(
    seq_len(nrow(data)),
    data$gene_entrez_id
  )
  output <- lapply(split_index, function(index) {
    group <- data[index, , drop = FALSE]
    constant_columns <- c(
      "gene_symbol", "gene_name", "gene_type",
      "gene_ensembl_id", "mrna_logFC",
      "mrna_P.Value", "mrna_adj.P.Val"
    )
    for (column in constant_columns) {
      values <- unique(as.character(group[[column]]))
      if (length(values) != 1L) {
        stop(
          column,
          " varies within a selected Entrez gene.",
          call. = FALSE
        )
      }
    }
    data.frame(
      gene_entrez_id = group$gene_entrez_id[1L],
      gene_symbol = group$gene_symbol[1L],
      gene_name = group$gene_name[1L],
      gene_type = group$gene_type[1L],
      gene_ensembl_id = group$gene_ensembl_id[1L],
      mrna_logFC = group$mrna_logFC[1L],
      mrna_P.Value = group$mrna_P.Value[1L],
      mrna_adj.P.Val = group$mrna_adj.P.Val[1L],
      selected_pair_count = nrow(group),
      selected_mirna_count = length(unique(group$mirna_id)),
      selected_mirna_ids = collapse_values(group$mirna_id),
      inverse_de_pair_count = sum(group$inverse_logFC_direction),
      compatible_unconfirmed_pair_count = sum(
        group$directionally_compatible_unconfirmed
      ),
      target_evidence_categories = collapse_values(
        group$target_evidence_category
      ),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, output)
  output <- output[
    order(output$mrna_adj.P.Val, output$gene_symbol),
    ,
    drop = FALSE
  ]
  row.names(output) <- NULL
  output
}

selected_genes <- collapse_selected_genes(integrated_pairs)
inverse_genes <- collapse_selected_genes(
  integrated_pairs[
    integrated_pairs$inverse_logFC_direction,
    ,
    drop = FALSE
  ]
)
compatible_genes <- collapse_selected_genes(compatible_pairs)

if (
  nrow(selected_genes) != 53L ||
    nrow(inverse_genes) != 12L ||
    nrow(compatible_genes) != 9L ||
    !all(
      selected_genes$gene_entrez_id %in%
        candidate_target_background$gene_entrez_id
    ) ||
    !all(
      inverse_genes$gene_entrez_id %in%
        selected_genes$gene_entrez_id
    ) ||
    !all(
      compatible_genes$gene_entrez_id %in%
        inverse_genes$gene_entrez_id
    )
) {
  stop(
    paste(
      "The formal ORA list or descriptive review lists failed",
      "expected size or nesting checks."
    ),
    call. = FALSE
  )
}

candidate_target_rows <- unique_pairs[
  unique_pairs$mirna_id %in% candidate_mirnas,
  ,
  drop = FALSE
]
target_link_key <- paste(
  candidate_target_rows$mirna_id,
  candidate_target_rows$gene_entrez_id,
  sep = "||"
)
target_link_groups <- split(
  seq_len(nrow(candidate_target_rows)),
  target_link_key
)
target_set_long <- lapply(target_link_groups, function(index) {
  group <- candidate_target_rows[index, , drop = FALSE]
  annotation_columns <- c(
    "mirna_id", "gene_entrez_id", "gene_symbol",
    "gene_name", "gene_type", "gene_ensembl_id"
  )
  for (column in annotation_columns) {
    if (length(unique(as.character(group[[column]]))) != 1L) {
      stop(
        column,
        " varies within a candidate miRNA--gene target set.",
        call. = FALSE
      )
    }
  }
  mirtarbase_supported <- any(group$mirtarbase_supported)
  targetscanmouse_supported <- any(
    group$targetscanmouse_supported
  )
  data.frame(
    mirna_id = group$mirna_id[1L],
    gene_entrez_id = group$gene_entrez_id[1L],
    gene_symbol = group$gene_symbol[1L],
    gene_name = group$gene_name[1L],
    gene_type = group$gene_type[1L],
    gene_ensembl_id = group$gene_ensembl_id[1L],
    mirtarbase_supported = mirtarbase_supported,
    targetscanmouse_supported = targetscanmouse_supported,
    target_evidence_category = if (
      mirtarbase_supported && targetscanmouse_supported
    ) {
      "miRTarBase_and_TargetScanMouse"
    } else if (mirtarbase_supported) {
      "miRTarBase_only"
    } else {
      "TargetScanMouse_only"
    },
    stringsAsFactors = FALSE
  )
})
target_set_long <- do.call(rbind, target_set_long)
target_set_long <- target_set_long[
  order(target_set_long$mirna_id, target_set_long$gene_entrez_id),
  ,
  drop = FALSE
]
row.names(target_set_long) <- NULL

target_set_sizes <- as.data.frame(
  table(target_set_long$mirna_id),
  stringsAsFactors = FALSE
)
names(target_set_sizes) <- c("mirna_id", "measured_target_gene_count")
target_set_sizes$measured_target_gene_count <- as.integer(
  target_set_sizes$measured_target_gene_count
)
target_set_sizes <- target_set_sizes[
  order(target_set_sizes$mirna_id),
  ,
  drop = FALSE
]
row.names(target_set_sizes) <- NULL

expected_target_set_sizes <- c(
  `mmu-mir-100-5p` = 24L,
  `mmu-mir-146a-5p` = 253L,
  `mmu-mir-1839-3p` = 139L,
  `mmu-mir-184-3p` = 29L,
  `mmu-mir-223-3p` = 557L,
  `mmu-mir-5130` = 157L
)
observed_target_set_sizes <- setNames(
  target_set_sizes$measured_target_gene_count,
  target_set_sizes$mirna_id
)
if (
  nrow(target_set_long) != 1159L ||
    anyDuplicated(paste(
      target_set_long$mirna_id,
      target_set_long$gene_entrez_id,
      sep = "||"
    )) ||
    !identical(
      observed_target_set_sizes[names(expected_target_set_sizes)],
      expected_target_set_sizes
    ) ||
    !all(
      target_set_long$gene_entrez_id %in%
        all_measured_background$gene_entrez_id
    )
) {
  stop(
    "Candidate-miRNA measured target gene sets failed validation.",
    call. = FALSE
  )
}

gmt_lines <- vapply(
  candidate_mirnas,
  function(mirna_id) {
    gene_ids <- sort(unique(
      target_set_long$gene_entrez_id[
        target_set_long$mirna_id == mirna_id
      ]
    ))
    paste(
      c(
        mirna_id,
        "measured_mouse_targets_miRTarBase_or_TargetScanMouse",
        gene_ids
      ),
      collapse = "\t"
    )
  },
  character(1)
)

dir.create(ora_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(review_list_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(background_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(target_set_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  all_measured_background,
  all_measured_background_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  all_target_linked_background,
  all_target_linked_background_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  candidate_target_background,
  candidate_target_background_file,
  row.names = FALSE,
  na = ""
)
write.csv(selected_genes, selected_gene_file, row.names = FALSE, na = "")
write.csv(inverse_genes, inverse_gene_file, row.names = FALSE, na = "")
write.csv(
  compatible_genes,
  compatible_gene_file,
  row.names = FALSE,
  na = ""
)
writeLines(
  c(
    "Descriptive small gene lists",
    "",
    paste(
      "The 12 inverse-DE genes and 9 directionally compatible but",
      "association-unconfirmed genes are retained for transparent",
      "review and annotation only."
    ),
    paste(
      "They are not formal ORA inputs because lists this small have",
      "low statistical power and can produce unstable results driven",
      "by one or two genes."
    ),
    paste(
      "Use full-ranked GSEA as the primary pathway analysis and the",
      "53-gene selected list as the only formal multi-omics ORA input."
    )
  ),
  review_note_file
)
write.csv(
  target_set_long,
  target_set_long_file,
  row.names = FALSE,
  na = ""
)
write.csv(
  target_set_sizes,
  target_set_size_file,
  row.names = FALSE,
  na = ""
)
writeLines(gmt_lines, target_set_gmt_file)

duplicated_entrez_gene_count <- sum(features_per_entrez > 1L)
redundant_transcript_cluster_count <- nrow(reference_de) -
  nrow(selected_features)
summary_lines <- c(
  "Multi-omics enrichment-input preparation summary",
  "",
  "This script prepares enrichment inputs and does not run enrichment tests.",
  "",
  paste0("mRNA DE model: treatment + sex + age"),
  paste0("mRNA contrasts prepared: ", length(contrast_order)),
  paste0(
    "mRNA transcript clusters before gene-level selection: ",
    nrow(reference_de)
  ),
  paste0(
    "Unique Entrez genes after gene-level selection: ",
    nrow(selected_features)
  ),
  paste0(
    "Entrez genes represented by multiple transcript clusters: ",
    duplicated_entrez_gene_count
  ),
  paste0(
    "Redundant transcript-cluster rows not selected: ",
    redundant_transcript_cluster_count
  ),
  paste(
    "Gene-level selection rule:",
    "highest AveExpr, then transcript-cluster ID for ties"
  ),
  "GSEA rank metric: limma moderated t-statistic",
  ""
)
for (contrast in contrast_order) {
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ", contrast, ": 20,082 ranked genes; exact rank-metric ties = ",
      gsea_tie_counts[[contrast]]
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  paste0(
    "All measured unique mRNA-gene background: ",
    nrow(all_measured_background)
  ),
  paste0(
    "All measured target-linked gene background: ",
    nrow(all_target_linked_background)
  ),
  paste0(
    paste(
      "Measured targets of the six selected miRNAs",
      "(formal multi-omics ORA background):"
    ),
    " ",
    nrow(candidate_target_background)
  ),
  "",
  paste0(
    "Formal ORA selected target-linked genes: ",
    nrow(selected_genes)
  ),
  paste0(
    "Descriptive inverse-DE review genes (not formal ORA): ",
    nrow(inverse_genes)
  ),
  paste0(
    paste(
      "Descriptive directionally compatible but unconfirmed",
      "review genes (not formal ORA):"
    ),
    " ",
    nrow(compatible_genes)
  ),
  paste(
    "Small-list decision: the 12- and 9-gene lists are retained only",
    "for descriptive annotation because standalone ORA would be",
    "underpowered and highly sensitive to one or two genes."
  ),
  "",
  "Candidate-miRNA measured target gene-set sizes:"
)
for (i in seq_len(nrow(target_set_sizes))) {
  summary_lines <- c(
    summary_lines,
    paste0(
      "- ", target_set_sizes$mirna_id[i], ": ",
      target_set_sizes$measured_target_gene_count[i]
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  paste0("Wrote GSEA files under: ", gsea_dir),
  paste0("Wrote formal ORA input under: ", ora_dir),
  paste0(
    "Wrote descriptive, non-ORA review lists under: ",
    review_list_dir
  ),
  paste0("Wrote backgrounds under: ", background_dir),
  paste0("Wrote miRNA target gene sets under: ", target_set_dir),
  "",
  paste(
    "Interpretation boundary: these files define tested gene lists,",
    "rankings, and backgrounds. They contain no pathway-enrichment",
    "p-values or biological conclusions."
  )
)
writeLines(summary_lines, summary_output_file)

cat(
  "GSEA rankings: ",
  length(contrast_order),
  " contrasts x ",
  nrow(selected_features),
  " unique Entrez genes\n",
  sep = ""
)
cat(
  "Formal multi-omics ORA background genes: ",
  nrow(candidate_target_background),
  "\n",
  sep = ""
)
cat(
  "Formal ORA genes: ",
  nrow(selected_genes),
  "; descriptive review genes: inverse=",
  nrow(inverse_genes),
  ", compatible=",
  nrow(compatible_genes),
  "\n",
  sep = ""
)
cat(
  "Candidate-miRNA target sets: ",
  nrow(target_set_sizes),
  " sets and ",
  nrow(target_set_long),
  " unique miRNA--gene links\n",
  sep = ""
)
cat("Wrote: ", summary_output_file, "\n", sep = "")
