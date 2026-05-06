#!/usr/bin/env Rscript

# Prepare one miRNA expression dataset.
# Usage:
#   Rscript scripts/mirna/02_prepare_expression_mirna.R vendor_chp
#   Rscript scripts/mirna/02_prepare_expression_mirna.R complete_rma_excluding_two
#   Rscript scripts/mirna/02_prepare_expression_mirna.R rma_normalized_mirna

args <- commandArgs(trailingOnly = TRUE)
valid_datasets <- c(
  "vendor_chp",
  "complete_rma_excluding_two",
  "rma_normalized_mirna"
)
if (length(args) != 1L || !args[[1]] %in% valid_datasets) {
  stop(
    "Supply exactly one dataset: ",
    paste(valid_datasets, collapse = ", "),
    call. = FALSE
  )
}
dataset <- args[[1]]

sample_sheet_file <- file.path(
  "results", "mirna", "sample_sheet", "sample_sheet_mirna.csv"
)
cel_dir <- file.path("data", "raw", "mirna_150001", "CEL Files")
chp_dir <- file.path("data", "processed", "mirna_150001", "cc-chp")
output_dir <- file.path("results", "mirna", "expression", dataset)
expression_file <- file.path(output_dir, "expression_matrix_mirna.csv")
expression_set_file <- file.path(output_dir, "expression_set_mirna.rds")
summary_file <- file.path(output_dir, "normalization_summary_mirna.txt")

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Missing required package: ", package, call. = FALSE)
  }
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

write_probe_matrix <- function(matrix_object, path, feature_name = NULL) {
  output <- data.frame(
    ProbeSetName = rownames(matrix_object),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!is.null(feature_name)) {
    output$feature_name <- feature_name
  }
  output <- cbind(output, as.data.frame(matrix_object, check.names = FALSE))
  write.csv(output, path, row.names = FALSE, quote = TRUE)
}

matrix_summary <- function(matrix_object) {
  values <- as.vector(matrix_object)
  c(
    paste("Dimensions:", paste(dim(matrix_object), collapse = " x ")),
    paste("Finite values:", sum(is.finite(values))),
    paste("NA values:", sum(is.na(values))),
    paste("NaN values:", sum(is.nan(values))),
    paste("Infinite values:", sum(is.infinite(values))),
    paste("Range:", paste(range(values, finite = TRUE), collapse = " to "))
  )
}

feature_names_from_eset <- function(eset) {
  ids <- Biobase::featureNames(eset)
  feature_names <- ids
  annotation_data <- tryCatch(
    oligo::getNetAffx(eset, "transcript"),
    error = function(error) NULL
  )
  if (!is.null(annotation_data)) {
    annotation_data <- as.data.frame(annotation_data)
    id_column <- intersect(c("transcriptclusterid", "probesetid"), names(annotation_data))
    name_column <- intersect(c("probesetid", "transcriptclusterid"), names(annotation_data))
    if (length(id_column) > 0L && length(name_column) > 0L) {
      matched <- match(ids, as.character(annotation_data[[id_column[[1]]]]))
      candidate <- as.character(annotation_data[[name_column[[1]]]][matched])
      feature_names[!is.na(candidate)] <- candidate[!is.na(candidate)]
    }
  }
  feature_names
}

prepare_vendor_chp <- function() {
  require_package("affxparser")
  stop_if_missing(sample_sheet_file, "Sample sheet")
  stop_if_missing(chp_dir, "Vendor CHP directory")

  sample_sheet <- read.csv(
    sample_sheet_file, stringsAsFactors = FALSE, check.names = FALSE
  )
  chp_files <- list.files(
    chp_dir,
    pattern = "[.]rma-dabg[.]chp$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(chp_files) == 0L) {
    stop("No vendor RMA-DABG CHP files found.", call. = FALSE)
  }

  cel_stem <- sub("[.]CEL$", "", sample_sheet$file_name, ignore.case = TRUE)
  chp_stem <- sub(
    "[.]rma-dabg[.]chp$", "", basename(chp_files), ignore.case = TRUE
  )
  match_index <- match(cel_stem, chp_stem)
  if (anyNA(match_index)) {
    stop("Not every sample-sheet CEL has an exactly matching CHP stem.", call. = FALSE)
  }

  keep <- sample_sheet$match_status == "matched" & !sample_sheet$is_control
  analysis_samples <- sample_sheet[keep, , drop = FALSE]
  analysis_files <- chp_files[match_index[keep]]

  read_entries <- function(path) {
    entries <- affxparser::readChp(path)$QuantificationEntries
    required <- c("ProbeSetName", "QuantificationValue", "PValue", "ID")
    if (is.null(entries) || !all(required %in% names(entries))) {
      stop("CHP lacks required quantification columns: ", path, call. = FALSE)
    }
    as.data.frame(entries[, required], stringsAsFactors = FALSE)
  }

  entries <- lapply(analysis_files, read_entries)
  reference_ids <- as.character(entries[[1]]$ProbeSetName)
  if (any(vapply(
    entries,
    function(x) !identical(as.character(x$ProbeSetName), reference_ids),
    logical(1)
  ))) {
    stop("CHP files do not have identical ProbeSetName order.", call. = FALSE)
  }

  sample_names <- cel_stem[keep]
  expression_matrix <- vapply(
    entries, function(x) as.numeric(x$QuantificationValue), numeric(length(reference_ids))
  )
  dabg_matrix <- vapply(
    entries, function(x) as.numeric(x$PValue), numeric(length(reference_ids))
  )
  rownames(expression_matrix) <- rownames(dabg_matrix) <- reference_ids
  colnames(expression_matrix) <- colnames(dabg_matrix) <- sample_names

  if (any(!is.finite(expression_matrix)) || any(!is.finite(dabg_matrix))) {
    stop("Vendor expression or DABG matrix contains non-finite values.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_probe_matrix(expression_matrix, expression_file)
  write_probe_matrix(
    dabg_matrix,
    file.path(output_dir, "dabg_pvalues_mirna.csv")
  )
  annotation <- unique(entries[[1]][, c("ProbeSetName", "ID")])
  write.csv(
    annotation,
    file.path(output_dir, "probe_annotation_mirna.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  write.csv(
    analysis_samples,
    file.path(output_dir, "analysis_samples_mirna.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  writeLines(
    c(
      "Dataset: vendor_chp",
      paste("CHP files used:", length(analysis_files)),
      matrix_summary(expression_matrix)
    ),
    summary_file
  )
}

prepare_cel_rma <- function(background) {
  require_package("Biobase")
  require_package("oligo")
  require_package("pd.mirna.4.1")
  stop_if_missing(cel_dir, "CEL directory")

  cel_files <- list.files(
    cel_dir, pattern = "[.]CEL$", full.names = TRUE, ignore.case = TRUE
  )
  if (length(cel_files) == 0L) {
    stop("No CEL files found.", call. = FALSE)
  }

  excluded <- character()
  if (background) {
    excluded <- c(
      "5505014461569112923888_A07_499422-46.CEL",
      "5505014461569112923889_G05_529725-4.CEL"
    )
    missing_exclusions <- setdiff(excluded, basename(cel_files))
    if (length(missing_exclusions) > 0L) {
      stop(
        "Expected problematic CEL file(s) not found: ",
        paste(missing_exclusions, collapse = ", "),
        call. = FALSE
      )
    }
    cel_files <- cel_files[!basename(cel_files) %in% excluded]
  }

  raw_cel <- oligo::read.celfiles(cel_files)
  normalized <- oligo::rma(raw_cel, background = background)
  expression_matrix <- Biobase::exprs(normalized)
  if (any(!is.finite(expression_matrix))) {
    stop("RMA produced non-finite expression values.", call. = FALSE)
  }

  colnames(expression_matrix) <- sub(
    "[.]CEL$", "", basename(colnames(expression_matrix)), ignore.case = TRUE
  )
  feature_names <- feature_names_from_eset(normalized)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_probe_matrix(expression_matrix, expression_file, feature_names)
  saveRDS(normalized, expression_set_file)
  if (length(excluded) > 0L) {
    write.csv(
      data.frame(file_name = excluded),
      file.path(output_dir, "excluded_arrays_mirna.csv"),
      row.names = FALSE,
      quote = TRUE
    )
  }
  writeLines(
    c(
      paste("Dataset:", dataset),
      paste("RMA background correction:", background),
      paste("CEL files used:", length(cel_files)),
      paste("CEL files excluded:", length(excluded)),
      if (length(excluded) > 0L) paste("Excluded:", paste(excluded, collapse = ", ")),
      matrix_summary(expression_matrix)
    ),
    summary_file
  )
}

if (dataset == "vendor_chp") {
  prepare_vendor_chp()
} else if (dataset == "complete_rma_excluding_two") {
  prepare_cel_rma(background = TRUE)
} else {
  prepare_cel_rma(background = FALSE)
}

message("Prepared dataset: ", dataset)
message("Expression matrix: ", expression_file)
