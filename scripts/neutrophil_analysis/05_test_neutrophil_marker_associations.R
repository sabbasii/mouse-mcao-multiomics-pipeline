# Test adjusted paired-animal associations for neutrophil-marker candidates
#
# Model: mRNA expression ~ miRNA expression + treatment + sex + age_group
# Associations are tested once per unique measured miRNA--marker pair across
# the 43 paired animals. This is supportive evidence only, not mediation or
# direct-regulation evidence.
#
# Inputs:
#   results/neutrophil_analysis/tables/neutrophil_marker_integrated_pairs.csv
#   results/multiomics/sample_manifest/paired_manifest_mirna_mrna.csv
#   results/mirna/expression/rma_normalized_mirna/annotation/mouse_mature_mirna_expression_mirna.csv
#   results/mrna/analysis_ready/expression_matrix_unique_gene_mapped_mrna.csv
#   results/mrna/analysis_ready/transcript_cluster_annotation_unique_gene_mapped_mrna.csv
#
# Outputs:
#   results/neutrophil_analysis/tables/neutrophil_marker_association_results.csv
#   results/neutrophil_analysis/tables/neutrophil_marker_association_summary.txt

options(stringsAsFactors = FALSE)
root <- normalizePath(".", mustWork = TRUE)
integrated_file <- file.path(root, "results/neutrophil_analysis/tables/neutrophil_marker_integrated_pairs.csv")
manifest_file <- file.path(root, "results/multiomics/sample_manifest/paired_manifest_mirna_mrna.csv")
mirna_file <- file.path(root, "results/mirna/expression/rma_normalized_mirna/annotation/mouse_mature_mirna_expression_mirna.csv")
mrna_file <- file.path(root, "results/mrna/analysis_ready/expression_matrix_unique_gene_mapped_mrna.csv")
mrna_annotation_file <- file.path(root, "results/mrna/analysis_ready/transcript_cluster_annotation_unique_gene_mapped_mrna.csv")
output_dir <- file.path(root, "results/neutrophil_analysis/tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required <- c(integrated_file, manifest_file, mirna_file, mrna_file, mrna_annotation_file)
if (any(!file.exists(required))) stop("Missing input: ", paste(required[!file.exists(required)], collapse = ", "), call. = FALSE)
read_table <- function(x) read.csv(x, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA", "N/A"))

candidates <- read_table(integrated_file)
manifest <- read_table(manifest_file)
mirna <- read_table(mirna_file)
mrna <- read_table(mrna_file)
annotation <- read_table(mrna_annotation_file)
if (nrow(manifest) != 43L) stop("Expected 43 paired animals; found ", nrow(manifest), ".", call. = FALSE)
needed_manifest <- c("animal_id", "treatment", "sex", "age_group", "mirna_file_name", "mrna_file_name")
if (length(setdiff(needed_manifest, names(manifest))) > 0) stop("Manifest lacks required columns.", call. = FALSE)
needed_candidates <- c("mirna_probeset_id", "mirna_id", "target_gene_symbol", "evidence_sources")
if (length(setdiff(needed_candidates, names(candidates))) > 0) stop("Integrated pairs lack required columns.", call. = FALSE)

annotation$SYMBOL <- trimws(as.character(annotation$SYMBOL))
annotation <- annotation[!duplicated(annotation$SYMBOL), c("SYMBOL", "transcript_cluster_id")]
candidates <- merge(candidates, annotation, by.x = "target_gene_symbol", by.y = "SYMBOL", all.x = FALSE, sort = FALSE)
candidates$pair_key <- paste(candidates$mirna_probeset_id, candidates$target_gene_symbol, sep = "||")
candidate_metadata <- aggregate(
  cbind(mirna_id, evidence_sources, contrast) ~ pair_key + mirna_probeset_id + target_gene_symbol + transcript_cluster_id,
  data = candidates,
  FUN = function(x) paste(sort(unique(as.character(x))), collapse = ";")
)

manifest$treatment <- factor(manifest$treatment, levels = c("Sham", "MCAO1hr", "MCAO3hr"))
manifest$sex <- factor(manifest$sex, levels = c("Female", "Male"))
manifest$age_group <- factor(manifest$age_group, levels = c("Young", "Old"))
if (!all(manifest$mirna_file_name %in% names(mirna)) || !all(manifest$mrna_file_name %in% names(mrna))) stop("Manifest samples are missing from an expression matrix.", call. = FALSE)
mirna$ProbeSetName <- as.character(mirna$ProbeSetName)
mrna$transcript_cluster_id <- as.character(mrna$transcript_cluster_id)

fit_one <- function(i) {
  m <- candidate_metadata[i, ]
  mirna_row <- match(m$mirna_probeset_id, mirna$ProbeSetName)
  mrna_row <- match(m$transcript_cluster_id, mrna$transcript_cluster_id)
  if (is.na(mirna_row) || is.na(mrna_row)) stop("Candidate feature is absent from an expression matrix.", call. = FALSE)
  d <- data.frame(
    mrna_expression = as.numeric(mrna[mrna_row, manifest$mrna_file_name]),
    mirna_expression = as.numeric(mirna[mirna_row, manifest$mirna_file_name]),
    treatment = manifest$treatment, sex = manifest$sex, age_group = manifest$age_group
  )
  if (any(!is.finite(as.matrix(d[, 1:2]))) || sd(d$mrna_expression) == 0 || sd(d$mirna_expression) == 0) stop("Non-finite or constant expression in candidate.", call. = FALSE)
  fit <- lm(mrna_expression ~ mirna_expression + treatment + sex + age_group, data = d)
  co <- summary(fit)$coefficients["mirna_expression", ]
  data.frame(
    n_paired_animals = nrow(d), unadjusted_pearson_r = cor(d$mirna_expression, d$mrna_expression),
    adjusted_beta = unname(co["Estimate"]), adjusted_SE = unname(co["Std. Error"]),
    adjusted_t = unname(co["t value"]), adjusted_P.Value = unname(co["Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
}
stats <- do.call(rbind, lapply(seq_len(nrow(candidate_metadata)), fit_one))
results <- cbind(candidate_metadata, stats)
results$adjusted_adj.P.Val <- p.adjust(results$adjusted_P.Value, method = "BH")
results$association_direction <- ifelse(results$adjusted_beta < 0, "negative", "positive")
results$adjusted_fdr_lt_0.10 <- results$adjusted_adj.P.Val < 0.10
results$adjusted_fdr_lt_0.05 <- results$adjusted_adj.P.Val < 0.05
results <- results[order(results$adjusted_P.Value), ]
rownames(results) <- NULL
write.csv(results, file.path(output_dir, "neutrophil_marker_association_results.csv"), row.names = FALSE, quote = TRUE)

summary_lines <- c(
  "Adjusted paired-animal neutrophil-marker association summary", "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Unique measured pairs tested:", nrow(results)),
  "Model: mRNA expression ~ miRNA expression + treatment + sex + age",
  paste("Negative adjusted associations:", sum(results$association_direction == "negative")),
  paste("Adjusted association p < 0.05:", sum(results$adjusted_P.Value < 0.05)),
  paste("Adjusted association FDR < 0.10:", sum(results$adjusted_fdr_lt_0.10)),
  paste("Adjusted association FDR < 0.05:", sum(results$adjusted_fdr_lt_0.05)),
  "",
  "Associations are exploratory and do not establish direct regulation,",
  "causality, or mediation. Target evidence remains database-derived."
)
writeLines(summary_lines, file.path(output_dir, "neutrophil_marker_association_summary.txt"))
message("Neutrophil marker association analysis complete. Pairs tested: ", nrow(results))
