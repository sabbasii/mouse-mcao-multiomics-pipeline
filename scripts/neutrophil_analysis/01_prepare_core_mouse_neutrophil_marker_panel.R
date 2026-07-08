# Prepare a source-traceable core mouse neutrophil marker panel
#
# This script freezes the marker panel before inspecting treatment-associated
# differential expression. It does not select genes because they are DE in the
# present experiment. The panel is based on independently profiled mouse
# neutrophils from ImmGen and sorted-neutrophil transcriptome studies.
#
# Inputs:
#   results/mrna/analysis_ready/transcript_cluster_annotation_unique_gene_mapped_mrna.csv
#     Measured, uniquely gene-mapped mRNA features used to assess coverage.
#   results/mrna/analysis_ready/transcript_cluster_annotation_unique_gene_mapped_mrna.csv
#     Measured mRNA gene universe used for the coverage check.
#
# Outputs:
#   results/neutrophil_analysis/gene_sets/core_mouse_neutrophil_marker_panel.csv
#     Frozen marker list with evidence tier and source provenance.
#   results/neutrophil_analysis/gene_sets/core_mouse_neutrophil_marker_coverage.csv
#     Measured/unmeasured status and identifiers for every marker.
#   results/neutrophil_analysis/gene_sets/core_mouse_neutrophil_marker_metadata.txt
#     Panel rationale, references, and coverage counts.

options(stringsAsFactors = FALSE)

project_root <- normalizePath(".", mustWork = TRUE)
annotation_path <- file.path(
  project_root,
  "results/mrna/analysis_ready/transcript_cluster_annotation_unique_gene_mapped_mrna.csv"
)
output_dir <- file.path(project_root, "results/neutrophil_analysis/gene_sets")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_file <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  TRUE
}
invisible(required_file(annotation_path))

annotation <- read.csv(annotation_path, check.names = FALSE)
needed_annotation <- c("ENTREZID", "SYMBOL", "GENENAME")
missing_annotation <- setdiff(needed_annotation, names(annotation))
if (length(missing_annotation) > 0) {
  stop("mRNA annotation is missing columns: ",
       paste(missing_annotation, collapse = ", "), call. = FALSE)
}

# Tier 1: lineage-associated genes repeatedly used to identify mouse
# neutrophils in independently profiled populations. Tier 2: genes that help
# represent maturation or inflammatory state and should not be interpreted as
# exclusive lineage markers.
panel <- data.frame(
  gene_symbol = c(
    "Ly6g", "S100a8", "S100a9", "Retnlg", "Ngp", "Lcn2", "Camp",
    "Mmp8", "Mmp9", "Cxcr2", "Fcgr3", "Itgam", "Csf3r", "Ltf",
    "Elane", "Mpo", "Il1b", "Cxcr4"
  ),
  panel_tier = c(
    rep("core_lineage", 16), rep("state_support", 2)
  ),
  evidence = c(
    "canonical Ly6G+ neutrophil identity marker",
    "mouse neutrophil-associated inflammatory marker",
    "mouse neutrophil-associated inflammatory marker",
    "blood/spleen neutrophil-associated marker",
    "mature-neutrophil-associated granule marker",
    "immature-neutrophil-associated marker; also inflammatory",
    "mouse neutrophil antimicrobial granule marker",
    "mouse neutrophil granule protease marker",
    "mouse neutrophil granule protease marker",
    "mouse neutrophil chemokine-receptor marker",
    "myeloid/neutrophil-associated Fc receptor marker",
    "CD11b myeloid adhesion marker; not neutrophil-exclusive",
    "neutrophil maturation receptor marker",
    "neutrophil granule marker",
    "neutrophil primary-granule protease marker",
    "neutrophil primary-granule peroxidase marker",
    "inflammatory neutrophil state marker",
    "aged/migratory neutrophil state marker"
  ),
  source_reference = c(
    rep("ImmGen neutrophil profiling; ImmGen single-cell neutrophils", 18)
  ),
  stringsAsFactors = FALSE
)

panel <- panel[!duplicated(panel$gene_symbol), ]
annotation$SYMBOL <- trimws(as.character(annotation$SYMBOL))
annotation$ENTREZID <- trimws(as.character(annotation$ENTREZID))
annotation$GENENAME <- trimws(as.character(annotation$GENENAME))
annotation <- annotation[!duplicated(annotation$SYMBOL), ]

coverage <- merge(
  panel,
  annotation[, c("SYMBOL", "ENTREZID", "GENENAME")],
  by.x = "gene_symbol", by.y = "SYMBOL", all.x = TRUE, sort = FALSE
)
coverage$measured_in_mrna_matrix <- !is.na(coverage$ENTREZID) & coverage$ENTREZID != ""
coverage <- coverage[order(match(coverage$gene_symbol, panel$gene_symbol)), ]
rownames(coverage) <- NULL

write.csv(panel, file.path(output_dir, "core_mouse_neutrophil_marker_panel.csv"),
          row.names = FALSE, quote = TRUE)
write.csv(coverage, file.path(output_dir, "core_mouse_neutrophil_marker_coverage.csv"),
          row.names = FALSE, quote = TRUE)

metadata_lines <- c(
  "Core mouse neutrophil marker panel",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Primary evidence sources:",
  "1. ImmGen-based mouse neutrophil profiling: PMC4184787.",
  "2. ImmGen single-cell neutrophil/neutrotime dataset: ImmGen Databrowser19.",
  "",
  "The core_lineage tier is intended for lineage-associated screening.",
  "The state_support tier represents maturation or inflammatory state and is not",
  "interpreted as neutrophil-exclusive.",
  "",
  paste("Markers in panel:", nrow(panel)),
  paste("Core-lineage markers:", sum(panel$panel_tier == "core_lineage")),
  paste("State-support markers:", sum(panel$panel_tier == "state_support")),
  paste("Measured in the analysis-ready mRNA annotation:", sum(coverage$measured_in_mrna_matrix)),
  paste("Not measured or not mapped:", sum(!coverage$measured_in_mrna_matrix)),
  "",
  "Bulk mRNA expression cannot distinguish neutrophil abundance from altered",
  "transcription within neutrophils. This panel is therefore a marker-based",
  "screen and not a cell-counting or cell-deconvolution result."
)
writeLines(metadata_lines, file.path(output_dir, "core_mouse_neutrophil_marker_metadata.txt"))

message("Core mouse neutrophil marker panel prepared.")
message("Markers: ", nrow(panel), "; measured: ", sum(coverage$measured_in_mrna_matrix))
