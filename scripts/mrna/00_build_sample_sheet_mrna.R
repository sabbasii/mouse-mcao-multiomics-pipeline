#!/usr/bin/env Rscript

# Build and validate the Clariom S/mRNA 150002 sample sheet.
# This script performs filename inventory and metadata matching only. It does
# not read CEL intensities, normalize arrays, or run differential expression.
#
# Run from the repository root:
#   Rscript scripts/mrna/00_build_sample_sheet_mrna.R

metadata_file <- file.path(
  "metadata",
  "Animal-Tracking-Sheet-RNA-Animal-ID.csv"
)
cel_dir <- file.path("data", "raw", "clariomS_150002", "CEL Files")
output_dir <- file.path("results", "mrna", "sample_sheet")

sample_sheet_file <- file.path(
  output_dir,
  "sample_sheet_mrna_150002.csv"
)
summary_file <- file.path(
  output_dir,
  "sample_matching_summary_mrna_150002.txt"
)
metadata_without_cel_file <- file.path(
  output_dir,
  "metadata_without_cel_mrna_150002.csv"
)

required_columns <- c(
  "Surgery Date",
  "Animal ID",
  "DOB",
  "Weight",
  "Sex",
  "Age",
  "Age Group",
  "Treatment"
)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

stop_if_missing(metadata_file, "Metadata file")

if (!dir.exists(cel_dir)) {
  stop("CEL directory not found: ", cel_dir, call. = FALSE)
}

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "N/A", "#N/A")
)

missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0L) {
  stop(
    "Metadata file is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

metadata$`Animal ID` <- trimws(metadata$`Animal ID`)

cel_files <- sort(list.files(
  cel_dir,
  pattern = "\\.CEL$",
  full.names = FALSE,
  ignore.case = TRUE
))

if (length(cel_files) == 0L) {
  stop("No CEL files found in: ", cel_dir, call. = FALSE)
}

if (anyDuplicated(tolower(cel_files))) {
  duplicated_files <- unique(cel_files[
    duplicated(tolower(cel_files)) |
      duplicated(tolower(cel_files), fromLast = TRUE)
  ])
  stop(
    "Duplicate CEL filenames found: ",
    paste(duplicated_files, collapse = ", "),
    call. = FALSE
  )
}

# Clariom S filenames begin with the animal ID followed by an underscore:
#   497727-35_5507404444215020923341_D07.CEL
# Control arrays begin with "Control_".
animal_id_from_filename <- sub("_.*$", "", cel_files)
is_control <- grepl("^Control_", cel_files, ignore.case = TRUE)
animal_id_from_filename[is_control] <- NA_character_

cel_animal_ids <- animal_id_from_filename[!is_control]

duplicate_metadata_ids <- unique(metadata$`Animal ID`[
  metadata$`Animal ID` %in% cel_animal_ids &
    (
      duplicated(metadata$`Animal ID`) |
        duplicated(metadata$`Animal ID`, fromLast = TRUE)
    )
])

if (length(duplicate_metadata_ids) > 0L) {
  stop(
    "Metadata contains duplicate rows for CEL animal ID(s): ",
    paste(duplicate_metadata_ids, collapse = ", "),
    ". Resolve these before building the sample sheet.",
    call. = FALSE
  )
}

cel_table <- data.frame(
  file_name = cel_files,
  cel_path = file.path(cel_dir, cel_files),
  animal_id = animal_id_from_filename,
  is_control = is_control,
  stringsAsFactors = FALSE
)

metadata_for_join <- metadata[, required_columns]
names(metadata_for_join) <- c(
  "surgery_date",
  "animal_id",
  "date_of_birth",
  "weight",
  "sex",
  "age",
  "age_group",
  "treatment"
)

# Blank metadata identifiers cannot match a biological CEL sample and would
# otherwise match the NA identifier used for control arrays by merge().
metadata_for_join <- metadata_for_join[
  !is.na(metadata_for_join$animal_id) &
    nzchar(metadata_for_join$animal_id),
  ,
  drop = FALSE
]

sample_sheet <- merge(
  cel_table,
  metadata_for_join,
  by = "animal_id",
  all.x = TRUE,
  sort = FALSE
)

# Restore the original CEL-file order after merge().
sample_sheet <- sample_sheet[
  match(cel_files, sample_sheet$file_name),
  ,
  drop = FALSE
]

sample_sheet$platform <- "clariomS_150002"
has_metadata_match <- !sample_sheet$is_control &
  sample_sheet$animal_id %in% metadata_for_join$animal_id
sample_sheet$match_status <- ifelse(
  sample_sheet$is_control,
  "control_without_animal_metadata",
  ifelse(has_metadata_match, "matched", "cel_without_metadata")
)

sample_sheet <- sample_sheet[, c(
  "platform",
  "file_name",
  "cel_path",
  "animal_id",
  "surgery_date",
  "date_of_birth",
  "weight",
  "sex",
  "age",
  "age_group",
  "treatment",
  "is_control",
  "match_status"
)]

metadata_without_cel <- metadata[
  !is.na(metadata$`Animal ID`) &
    nzchar(metadata$`Animal ID`) &
    !metadata$`Animal ID` %in% cel_animal_ids,
  required_columns,
  drop = FALSE
]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(sample_sheet, sample_sheet_file, row.names = FALSE, na = "")
write.csv(
  metadata_without_cel,
  metadata_without_cel_file,
  row.names = FALSE,
  na = ""
)

summary_lines <- c(
  "Clariom S/mRNA 150002 sample matching summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("Metadata file:", metadata_file),
  paste("CEL directory:", cel_dir),
  "",
  paste("Metadata rows:", nrow(metadata)),
  paste("CEL files:", nrow(cel_table)),
  paste("Control CEL files:", sum(sample_sheet$is_control)),
  paste("Matched non-control CEL files:", sum(
    sample_sheet$match_status == "matched"
  )),
  paste("Non-control CEL files without metadata:", sum(
    sample_sheet$match_status == "cel_without_metadata"
  )),
  paste("Metadata rows without a CEL file:", nrow(metadata_without_cel)),
  paste("Duplicate CEL filenames:", sum(duplicated(tolower(cel_files)))),
  paste("Duplicate matched metadata animal IDs:", length(
    duplicate_metadata_ids
  )),
  "",
  "Match-status counts:",
  capture.output(print(table(
    sample_sheet$match_status,
    useNA = "ifany"
  ))),
  "",
  "Sex counts for matched non-control samples:",
  capture.output(print(table(
    sample_sheet$sex[sample_sheet$match_status == "matched"],
    useNA = "ifany"
  ))),
  "",
  "Age-group counts for matched non-control samples:",
  capture.output(print(table(
    sample_sheet$age_group[sample_sheet$match_status == "matched"],
    useNA = "ifany"
  ))),
  "",
  "Treatment counts for matched non-control samples:",
  capture.output(print(table(
    sample_sheet$treatment[sample_sheet$match_status == "matched"],
    useNA = "ifany"
  ))),
  "",
  paste("Sample sheet output:", sample_sheet_file),
  paste("Metadata-without-CEL output:", metadata_without_cel_file),
  paste("Summary output:", summary_file)
)

writeLines(summary_lines, summary_file)

message("Wrote sample sheet: ", sample_sheet_file)
message("Wrote metadata-without-CEL report: ", metadata_without_cel_file)
message("Wrote matching summary: ", summary_file)
