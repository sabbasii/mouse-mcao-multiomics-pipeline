#!/usr/bin/env Rscript

# Build and validate the animal-level miRNA-mRNA paired sample manifest.
#
# Inputs:
#   - miRNA sample sheet for the selected rma_normalized_mirna dataset
#   - mRNA analysis-ready sample metadata
#   - headers of the corresponding analysis-ready expression matrices
#
# Outputs:
#   results/multiomics/sample_manifest/
#     paired_manifest_mirna_mrna.csv
#       Exact pairs retained for the Sham, MCAO1hr, and MCAO3hr analysis.
#     unpaired_samples_mirna_mrna.csv
#       Biological samples present in only one molecular layer.
#     paired_manifest_summary.txt
#       Input checks, group counts, exclusions, and validation results.
#
# Usage:
#   Rscript scripts/multiomics/00_build_paired_sample_manifest.R

mirna_sample_file <- file.path(
  "results", "mirna", "sample_sheet", "sample_sheet_mirna.csv"
)
mirna_expression_file <- file.path(
  "results", "mirna", "expression", "rma_normalized_mirna", "annotation",
  "mouse_mature_mirna_expression_mirna.csv"
)
mrna_sample_file <- file.path(
  "results", "mrna", "analysis_ready", "analysis_samples_mrna.csv"
)
mrna_expression_file <- file.path(
  "results", "mrna", "analysis_ready",
  "expression_matrix_unique_gene_mapped_mrna.csv"
)
output_dir <- file.path("results", "multiomics", "sample_manifest")

paired_file <- file.path(
  output_dir, "paired_manifest_mirna_mrna.csv"
)
unpaired_file <- file.path(
  output_dir, "unpaired_samples_mirna_mrna.csv"
)
summary_file <- file.path(
  output_dir, "paired_manifest_summary.txt"
)

required_inputs <- c(
  mirna_sample_file,
  mirna_expression_file,
  mrna_sample_file,
  mrna_expression_file
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    call. = FALSE
  )
}

mirna_samples <- read.csv(
  mirna_sample_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)
mrna_samples <- read.csv(
  mrna_sample_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A")
)

required_sample_columns <- c(
  "file_name", "animal_id", "sex", "age", "treatment",
  "is_control", "match_status"
)
for (layer in c("mirna", "mrna")) {
  sample_table <- get(paste0(layer, "_samples"))
  missing_columns <- setdiff(required_sample_columns, names(sample_table))
  if (length(missing_columns) > 0L) {
    stop(
      layer,
      " sample table lacks required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

mirna_samples <- mirna_samples[
  !mirna_samples$is_control &
    mirna_samples$match_status == "matched",
  ,
  drop = FALSE
]
mrna_samples <- mrna_samples[
  !mrna_samples$is_control &
    mrna_samples$match_status == "matched",
  ,
  drop = FALSE
]

mirna_samples$manifest_age_group <- mirna_samples$age
mrna_samples$manifest_age_group <- if (
  "age_group" %in% names(mrna_samples)
) {
  mrna_samples$age_group
} else {
  mrna_samples$age
}

if (
  anyNA(mirna_samples$animal_id) ||
    anyNA(mrna_samples$animal_id) ||
    anyDuplicated(mirna_samples$animal_id) ||
    anyDuplicated(mrna_samples$animal_id)
) {
  stop(
    "Biological sample tables contain a missing or duplicated animal ID.",
    call. = FALSE
  )
}

mirna_expression_columns <- names(read.csv(
  mirna_expression_file,
  nrows = 0L,
  check.names = FALSE
))
mrna_expression_columns <- names(read.csv(
  mrna_expression_file,
  nrows = 0L,
  check.names = FALSE
))

if (!all(mirna_samples$file_name %in% mirna_expression_columns)) {
  stop(
    "Not every matched miRNA sample is present in the selected expression matrix.",
    call. = FALSE
  )
}
if (!all(mrna_samples$file_name %in% mrna_expression_columns)) {
  stop(
    "Not every matched mRNA sample is present in the analysis-ready expression matrix.",
    call. = FALSE
  )
}

mirna_treatment_map <- c(
  Sham = "Sham",
  MCAO1hr = "MCAO1hr",
  MCAO3hr = "MCAO3hr",
  MCAO24hr = "MCAO24hr"
)
mrna_treatment_map <- c(
  Sham = "Sham",
  `MCAO1hr` = "MCAO1hr",
  `MCAO3hr` = "MCAO3hr",
  MCAO24hr = "MCAO24hr",
  Naive = "Naive"
)

standardize_treatment <- function(values, mapping, layer) {
  standardized <- unname(mapping[values])
  if (anyNA(standardized)) {
    stop(
      "Unrecognized ",
      layer,
      " treatment label(s): ",
      paste(unique(values[is.na(standardized)]), collapse = ", "),
      call. = FALSE
    )
  }
  standardized
}

mirna_samples$standardized_treatment <- standardize_treatment(
  mirna_samples$treatment,
  mirna_treatment_map,
  "miRNA"
)
mrna_samples$standardized_treatment <- standardize_treatment(
  mrna_samples$treatment,
  mrna_treatment_map,
  "mRNA"
)

paired_ids <- sort(intersect(
  mirna_samples$animal_id,
  mrna_samples$animal_id
))

mirna_paired <- mirna_samples[
  match(paired_ids, mirna_samples$animal_id),
  ,
  drop = FALSE
]
mrna_paired <- mrna_samples[
  match(paired_ids, mrna_samples$animal_id),
  ,
  drop = FALSE
]

if (!identical(mirna_paired$animal_id, mrna_paired$animal_id)) {
  stop("Animal order differs between paired sample tables.", call. = FALSE)
}
if (!identical(
  mirna_paired$standardized_treatment,
  mrna_paired$standardized_treatment
)) {
  disagreement <- paired_ids[
    mirna_paired$standardized_treatment !=
      mrna_paired$standardized_treatment
  ]
  stop(
    "Treatment disagreement between molecular layers for: ",
    paste(disagreement, collapse = ", "),
    call. = FALSE
  )
}
if (!identical(mirna_paired$sex, mrna_paired$sex)) {
  disagreement <- paired_ids[mirna_paired$sex != mrna_paired$sex]
  stop(
    "Sex disagreement between molecular layers for: ",
    paste(disagreement, collapse = ", "),
    call. = FALSE
  )
}

analysis_treatments <- c("Sham", "MCAO1hr", "MCAO3hr")

all_exact_pairs <- data.frame(
  animal_id = paired_ids,
  treatment = mirna_paired$standardized_treatment,
  sex = mirna_paired$sex,
  age_group = mirna_paired$age,
  mirna_file_name = mirna_paired$file_name,
  mrna_file_name = mrna_paired$file_name,
  mirna_treatment_original = mirna_paired$treatment,
  mrna_treatment_original = mrna_paired$treatment,
  mrna_surgery_date = if ("surgery_date" %in% names(mrna_paired)) {
    mrna_paired$surgery_date
  } else {
    NA_character_
  },
  stringsAsFactors = FALSE,
  check.names = FALSE
)

paired_manifest <- all_exact_pairs[
  all_exact_pairs$treatment %in% analysis_treatments,
  ,
  drop = FALSE
]

make_unpaired_rows <- function(sample_table, layer, other_ids) {
  unpaired <- sample_table[
    !sample_table$animal_id %in% other_ids,
    ,
    drop = FALSE
  ]
  if (nrow(unpaired) == 0L) {
    return(data.frame(
      animal_id = character(),
      available_layer = character(),
      treatment_original = character(),
      standardized_treatment = character(),
      sex = character(),
      age_group = character(),
      file_name = character(),
      reason = character(),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    animal_id = unpaired$animal_id,
    available_layer = layer,
    treatment_original = unpaired$treatment,
    standardized_treatment = unpaired$standardized_treatment,
    sex = unpaired$sex,
    age_group = unpaired$manifest_age_group,
    file_name = unpaired$file_name,
    reason = paste0("missing_", if (layer == "miRNA") "mRNA" else "miRNA"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

unpaired_samples <- rbind(
  make_unpaired_rows(mirna_samples, "miRNA", mrna_samples$animal_id),
  make_unpaired_rows(mrna_samples, "mRNA", mirna_samples$animal_id)
)
unpaired_samples <- unpaired_samples[
  order(unpaired_samples$available_layer, unpaired_samples$animal_id),
  ,
  drop = FALSE
]

expected_all_counts <- c(
  Sham = 16L,
  MCAO1hr = 16L,
  MCAO3hr = 11L,
  MCAO24hr = 1L
)
observed_all_counts <- table(factor(
  all_exact_pairs$treatment,
  levels = names(expected_all_counts)
))
expected_analysis_counts <- expected_all_counts[analysis_treatments]
observed_analysis_counts <- table(factor(
  paired_manifest$treatment,
  levels = names(expected_analysis_counts)
))

stopifnot(
  nrow(all_exact_pairs) == 44L,
  nrow(paired_manifest) == 43L,
  !anyDuplicated(paired_manifest$animal_id),
  !anyDuplicated(paired_manifest$mirna_file_name),
  !anyDuplicated(paired_manifest$mrna_file_name),
  identical(
    unname(as.integer(observed_all_counts)),
    unname(as.integer(expected_all_counts))
  ),
  identical(
    unname(as.integer(observed_analysis_counts)),
    unname(as.integer(expected_analysis_counts))
  ),
  all(paired_manifest$mirna_file_name %in% mirna_expression_columns),
  all(paired_manifest$mrna_file_name %in% mrna_expression_columns),
  !any(paired_manifest$treatment == "MCAO24hr"),
  sum(all_exact_pairs$treatment == "MCAO24hr") == 1L
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  paired_manifest,
  paired_file,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)
write.csv(
  unpaired_samples,
  unpaired_file,
  row.names = FALSE,
  quote = TRUE,
  na = ""
)

format_counts <- function(counts) {
  paste(
    paste0(names(counts), ": ", as.integer(counts)),
    collapse = ", "
  )
}

input_md5 <- tools::md5sum(required_inputs)

writeLines(
  c(
    "Paired miRNA-mRNA sample manifest",
    "",
    paste("miRNA sample input:", mirna_sample_file),
    paste("miRNA expression input:", mirna_expression_file),
    paste("mRNA sample input:", mrna_sample_file),
    paste("mRNA expression input:", mrna_expression_file),
    "",
    paste("Matched biological miRNA samples:", nrow(mirna_samples)),
    paste("Matched biological mRNA samples:", nrow(mrna_samples)),
    paste("Exact animal-level pairs before treatment filtering:", nrow(all_exact_pairs)),
    paste("Paired analysis samples:", nrow(paired_manifest)),
    paste("Exact pairs excluded by treatment:", nrow(all_exact_pairs) - nrow(paired_manifest)),
    paste("Unpaired biological samples:", nrow(unpaired_samples)),
    "",
    paste("All paired treatment counts:", format_counts(observed_all_counts)),
    paste(
      "Paired analysis treatment counts:",
      format_counts(observed_analysis_counts)
    ),
    "",
    "Analysis treatments: Sham, MCAO1hr, MCAO3hr",
    "The single MCAO24hr exact pair is excluded from the analysis manifest.",
    "Naive is mRNA-only and is absent from the paired manifest.",
    "",
    "Validation:",
    "- Animal IDs are unique within each layer and manifest.",
    "- Sex agrees between layers for every paired animal.",
    "- Standardized treatment agrees between layers for every paired animal.",
    "- Every paired filename exists in its selected expression matrix.",
    "- The paired analysis manifest contains 43 mice.",
    "",
    "Input MD5 checksums:",
    paste(names(input_md5), unname(input_md5), sep = ": ")
  ),
  summary_file
)

message("Wrote paired manifest: ", paired_file)
message("Wrote unpaired sample report: ", unpaired_file)
message("Wrote manifest summary: ", summary_file)
