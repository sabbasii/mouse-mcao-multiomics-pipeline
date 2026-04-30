#!/usr/bin/env Rscript

# Build and validate the miRNA 150001 sample sheet.
# This script performs sample matching only. It does not normalize arrays or run
# differential expression.

metadata_file <- file.path(
  "metadata",
  "Animal-Tracking-Sheet-RNA-R-miRNA.csv"
)
cel_dir <- file.path("data", "raw", "mirna_150001", "CEL Files")
output_dir <- file.path("results", "sample_sheet")
sample_sheet_file <- file.path(output_dir, "sample_sheet_mirna_150001.csv")
summary_file <- file.path(output_dir, "sample_matching_summary_mirna_150001.txt")

required_columns <- c("File Name", "Animal ID", "Sex", "Age", "Treatment")

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

stop_if_missing(metadata_file, "Metadata file")
stop_if_missing(cel_dir, "CEL directory")

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

missing_columns <- setdiff(required_columns, names(metadata))
if (length(missing_columns) > 0) {
  stop(
    "Metadata file is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

cel_files <- list.files(
  cel_dir,
  pattern = "\\.CEL$",
  full.names = FALSE,
  ignore.case = TRUE
)

cel_table <- data.frame(
  "File Name" = sort(cel_files),
  cel_path = file.path(cel_dir, sort(cel_files)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sample_sheet <- merge(
  cel_table,
  metadata[, required_columns],
  by = "File Name",
  all.x = TRUE,
  all.y = TRUE,
  sort = TRUE
)

sample_sheet$platform <- "mirna_150001"
sample_sheet$is_control <- sample_sheet$Sex == "N/A" |
  sample_sheet$Age == "N/A" |
  sample_sheet$Treatment == "N/A" |
  grepl("Control", sample_sheet$`File Name`, ignore.case = TRUE)

sample_sheet$match_status <- ifelse(
  is.na(sample_sheet$cel_path),
  "metadata_without_cel",
  ifelse(is.na(sample_sheet$`Animal ID`), "cel_without_metadata", "matched")
)

sample_sheet <- sample_sheet[
  ,
  c(
    "platform",
    "File Name",
    "cel_path",
    "Animal ID",
    "Sex",
    "Age",
    "Treatment",
    "is_control",
    "match_status"
  )
]

names(sample_sheet)[names(sample_sheet) == "File Name"] <- "file_name"
names(sample_sheet)[names(sample_sheet) == "Animal ID"] <- "animal_id"
names(sample_sheet)[names(sample_sheet) == "Sex"] <- "sex"
names(sample_sheet)[names(sample_sheet) == "Age"] <- "age"
names(sample_sheet)[names(sample_sheet) == "Treatment"] <- "treatment"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(sample_sheet, sample_sheet_file, row.names = FALSE)

summary_lines <- c(
  "miRNA 150001 sample matching summary",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  paste("Metadata file:", metadata_file),
  paste("CEL directory:", cel_dir),
  "",
  paste("Metadata rows:", nrow(metadata)),
  paste("CEL files:", nrow(cel_table)),
  paste("Matched rows:", sum(sample_sheet$match_status == "matched")),
  paste(
    "Metadata rows without CEL:",
    sum(sample_sheet$match_status == "metadata_without_cel")
  ),
  paste(
    "CEL files without metadata:",
    sum(sample_sheet$match_status == "cel_without_metadata")
  ),
  paste("Control/N/A rows flagged:", sum(sample_sheet$is_control, na.rm = TRUE)),
  "",
  "Sex counts:",
  capture.output(print(table(sample_sheet$sex, useNA = "ifany"))),
  "",
  "Age counts:",
  capture.output(print(table(sample_sheet$age, useNA = "ifany"))),
  "",
  "Treatment counts:",
  capture.output(print(table(sample_sheet$treatment, useNA = "ifany"))),
  "",
  paste("Sample sheet output:", sample_sheet_file),
  paste("Summary output:", summary_file)
)

writeLines(summary_lines, summary_file)

message("Wrote sample sheet: ", sample_sheet_file)
message("Wrote matching summary: ", summary_file)
