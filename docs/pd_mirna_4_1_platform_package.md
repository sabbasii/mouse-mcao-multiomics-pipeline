# Build a `pd.mirna.4.1` Platform Package for Affymetrix miRNA 4.1 Arrays

## Table of contents

- [1. Purpose and scope](#1-purpose-and-scope)
- [2. What you need before starting](#2-what-you-need-before-starting)
  - [2.1. A compatible R and Bioconductor environment](#21-a-compatible-r-and-bioconductor-environment)
  - [2.2. Manufacturer array-design files](#22-manufacturer-array-design-files)
  - [2.3. Biological annotation](#23-biological-annotation)
  - [2.4. Raw CEL files](#24-raw-cel-files)
  - [2.5. Local repository layout](#25-local-repository-layout)
- [3. Build workflow](#3-build-workflow)
- [4. Reproduce the project environment](#4-reproduce-the-project-environment)
- [5. Prepare and validate the external files](#5-prepare-and-validate-the-external-files)
  - [5.1. Define the input paths](#51-define-the-input-paths)
  - [5.2. Confirm the platform identifiers](#52-confirm-the-platform-identifiers)
  - [5.3. Confirm the geometry and design counts](#53-confirm-the-geometry-and-design-counts)
- [6. Build the platform package](#6-build-the-platform-package)
  - [6.1. Define package metadata](#61-define-package-metadata)
  - [6.2. Create the miRNA platform seed](#62-create-the-mirna-platform-seed)
  - [6.3. Generate the package source](#63-generate-the-package-source)
- [7. Validate the generated package database](#7-validate-the-generated-package-database)
- [8. Install and verify the package](#8-install-and-verify-the-package)
- [9. Validate the package with CEL files](#9-validate-the-package-with-cel-files)
- [10. Import the biological annotation](#10-import-the-biological-annotation)
- [11. Expected outputs and checkpoints](#11-expected-outputs-and-checkpoints)
- [12. Troubleshooting](#12-troubleshooting)
  - [12.1. External-file problems](#121-external-file-problems)
  - [12.2. Package-build problems](#122-package-build-problems)
  - [12.3. R environment and installation problems](#123-r-environment-and-installation-problems)
  - [12.4. CEL-import problems](#124-cel-import-problems)
  - [12.5. RMA-normalization problems](#125-rma-normalization-problems)
- [13. References](#13-references)

---

## 1. Purpose and scope

Affymetrix GeneChip miRNA 4.1 CEL files contain probe-level fluorescence
intensities. They do not contain the complete array design needed to connect
physical probe coordinates to probesets. When
[`oligo`](https://bioconductor.org/packages/release/bioc/html/oligo.html)
reads these CEL files, it expects an R platform-design package named
`pd.mirna.4.1`.

At the time this workflow was developed, a compatible package was not
available from the active Bioconductor repository. This guide builds it from
the official Probe Group File (PGF) and Chip Layout File (CLF) using
[`pdInfoBuilder`](https://bioconductor.org/packages/release/bioc/html/pdInfoBuilder.html).

The completed package provides:

- physical array coordinates and probe IDs;
- probeset membership;
- perfect-match probe records and sequences; and
- the mappings used by `oligo` to summarize CEL intensities.

The package does not contain biological miRNA names. That information comes
from a separate NetAffx annotation CSV and is added only after probeset IDs
exist.

This document therefore follows one dependency chain:

```text
R/Bioconductor environment
        +
official PGF and CLF
        ↓
generated pd.mirna.4.1 source and SQLite database
        ↓
installed pd.mirna.4.1 package
        +
raw CEL files
        ↓
validated probe-level data
        +
NetAffx annotation CSV
        ↓
biologically annotated probesets
```

---

## 2. What you need before starting

Gather the software and external files in this section before beginning the
build. Later sections introduce each item again at the point where it is used.

### 2.1. A compatible R and Bioconductor environment

The package is built and tested inside an R project environment. This
repository uses [`renv`](https://rstudio.github.io/renv/) to keep package
versions reproducible.

The relevant project files are:

- `.Rprofile` — activates `renv` when R starts in the repository;
- `renv/activate.R` — bootstraps the project library;
- `renv.lock` — records package sources and versions;
- `DESCRIPTION` — lists the packages required by the analysis.

When `renv` is active, packages are installed into a project-specific library
instead of relying only on the user’s global R library. This prevents an
analysis from silently using unrelated package versions installed for another
project.

Required packages include:

- [`BiocManager`](https://cran.r-project.org/package=BiocManager);
- [`pdInfoBuilder`](https://bioconductor.org/packages/release/bioc/html/pdInfoBuilder.html);
- [`oligo`](https://bioconductor.org/packages/release/bioc/html/oligo.html);
- `affxparser`;
- `oligoClasses`;
- `Biobase`;
- `DBI`;
- `RSQLite`;
- `methods`;
- `renv` when reproducing this repository.

R and Bioconductor releases must be compatible. Avoid installing the custom
package with one R version and running it with another.

### 2.2. Manufacturer array-design files

The two files required to build `pd.mirna.4.1` are:

- **PGF:** `miRNA-4_1-st-v1.pgf`
  - defines probesets and their member probes;
- **CLF:** `miRNA-4_1.clf`
  - maps probe IDs to physical `x` and `y` array coordinates.

Download **MiRNA 4.1 Analysis Files for Expression Console and AGCC** from
Thermo Fisher’s
[GeneChip Array Library Files](https://www.thermofisher.com/br/en/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-library-files.html)
page.

Use the Expression Console/AGCC analysis bundle. The TAC configuration file
is not a replacement for the PGF and CLF. This distinction is also noted in
the
[Bioconductor miRNA 4.1 support discussion](https://support.bioconductor.org/p/96882/).

If the legacy bundle is no longer listed, request the **miRNA 4.1 Expression
Console/AGCC analysis library files** from Thermo Fisher technical support.

### 2.3. Biological annotation

The downstream annotation file is:

- `miRNA-4_1-st-v1.annotations.20160922.csv`

Obtain it from Thermo Fisher’s miRNA 4.1 NetAffx annotation resources. It maps
`Probe Set ID` to miRNA names, accessions, sequence types, species, and other
fields.

This CSV is not used to build the platform package. It is introduced later,
after the package has produced probeset identifiers.

### 2.4. Raw CEL files

CEL files are experimental inputs supplied by the laboratory, institutional
archive, or public repository associated with a study. They are not included
in this GitHub repository or in Thermo Fisher’s platform-design bundle.

Only a representative CEL file is needed for an initial package test. The
complete experiment can be loaded after that test succeeds.

### 2.5. Local repository layout

Manufacturer files, CEL files, build artifacts, and analysis results are
local-only and excluded from Git. After the external files are supplied, the
relevant layout is:

```text
mirna-microarray-pipeline/
├── .Rprofile
├── DESCRIPTION
├── renv.lock
├── renv/
│   └── activate.R
├── data/
│   └── raw/
│       └── mirna_150001/
│           └── CEL Files/
│               └── *.CEL
├── resources/
│   └── affymetrix_library_files/
│       └── mirna_4_1/
│           ├── miRNA-4_1-st-v1.pgf
│           ├── miRNA-4_1.clf
│           ├── miRNA-4_1-st-v1.annotations.20160922.csv
│           └── miRNA_NetAffx-CSV-Files.README.txt
├── build/
│   └── pd_mirna_4_1/
├── docs/
│   └── pd_mirna_4_1_platform_package.md
└── scripts/
    └── mirna/
```

Users cloning the repository must create the local input directories and
supply files they are licensed or authorized to use.

---

## 3. Build workflow

The steps below are performed in order. Each step consumes an input prepared
earlier and produces something required by the next step.

```text
1. Restore and inspect the project R environment
↓
2. Place the PGF and CLF in the local resource directory
↓
3. Validate platform identifiers, geometry, and design counts
↓
4. Combine the validated files with package metadata in a seed object
↓
5. Run pdInfoBuilder to generate the R package source and SQLite database
↓
6. Validate the generated database against the PGF checkpoints
↓
7. Install the generated package into the active project library
↓
8. Read one CEL file, then the full experiment, with oligo
↓
9. Import the NetAffx CSV and join biological annotation by probeset ID
```

---

## 4. Reproduce the project environment

This first operational step ensures that subsequent commands use one
consistent R installation and project library.

Start R from the repository root. Confirm the working directory and active
library:

```r
getwd()
R.version.string
.libPaths()
```

If `renv` is available, inspect the project:

```r
if (requireNamespace("renv", quietly = TRUE)) {
  renv::project()
  renv::status()
}
```

Restore packages recorded by the lockfile:

```r
if (requireNamespace("renv", quietly = TRUE)) {
  renv::restore()
}
```

The restore may report that the locally generated `pd.mirna.4.1` package
cannot yet be restored from a public repository. That is expected at this
stage: this guide builds and installs it in Sections 6–8.

Install any missing build dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioconductor_packages <- c(
  "pdInfoBuilder",
  "oligo",
  "affxparser",
  "oligoClasses",
  "Biobase",
  "RSQLite"
)

missing_bioc <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc)
}

required_packages <- c(
  bioconductor_packages,
  "DBI",
  "methods"
)

available <- vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1)
)

available
stopifnot(all(available))
```

At the end of this section, the build tools are available in the same R
environment that will later install and test the generated package.

---

## 5. Prepare and validate the external files

This section starts with the downloaded PGF and CLF. It produces validated
file paths and expected design counts, which are used when creating and
checking the package.

### 5.1. Define the input paths

Create the local resource directory:

```bash
mkdir -p resources/affymetrix_library_files/mirna_4_1
```

Place the extracted manufacturer files there, then define their paths:

```r
library_dir <- file.path(
  "resources",
  "affymetrix_library_files",
  "mirna_4_1"
)

pgf_file <- file.path(library_dir, "miRNA-4_1-st-v1.pgf")
clf_file <- file.path(library_dir, "miRNA-4_1.clf")
annotation_file <- file.path(
  library_dir,
  "miRNA-4_1-st-v1.annotations.20160922.csv"
)

design_files <- c(pgf_file, clf_file)

stopifnot(all(file.exists(design_files)))
stopifnot(all(file.info(design_files)$size > 0))

file.info(design_files)[, c("size", "mtime"), drop = FALSE]
```

The output of this step is a pair of existing, non-empty design-file paths.

### 5.2. Confirm the platform identifiers

Read the manufacturer headers from the files prepared in Section 5.1:

```r
pgf_header <- readLines(pgf_file, n = 15)
clf_header <- readLines(clf_file, n = 15)

header_value <- function(lines, key) {
  match <- grep(paste0("^#%", key, "="), lines, value = TRUE)

  if (length(match) != 1L) {
    stop("Expected exactly one header value for: ", key)
  }

  sub(paste0("^#%", key, "="), "", match)
}

pgf_chip_type <- header_value(pgf_header, "chip_type")
clf_chip_type <- header_value(clf_header, "chip_type")

c(PGF = pgf_chip_type, CLF = clf_chip_type)

stopifnot(pgf_chip_type == "miRNA-4_1")
stopifnot(clf_chip_type == "miRNA-4_1")
```

Do not continue if either file identifies another platform.

### 5.3. Confirm the geometry and design counts

The platform identifiers now match. Next, confirm that the two files also
describe the same physical design:

```r
pgf_rows <- as.integer(header_value(pgf_header, "num-rows"))
pgf_columns <- as.integer(header_value(pgf_header, "num-cols"))
pgf_probesets <- as.integer(header_value(pgf_header, "probesets"))

clf_rows <- as.integer(header_value(clf_header, "rows"))
clf_columns <- as.integer(header_value(clf_header, "cols"))
clf_features <- as.integer(header_value(clf_header, "datalines"))

design_summary <- c(
  PGF_rows = pgf_rows,
  PGF_columns = pgf_columns,
  PGF_probesets = pgf_probesets,
  CLF_rows = clf_rows,
  CLF_columns = clf_columns,
  CLF_features = clf_features
)

design_summary

stopifnot(pgf_rows == 541L)
stopifnot(pgf_columns == 541L)
stopifnot(clf_rows == 541L)
stopifnot(clf_columns == 541L)
stopifnot(clf_rows * clf_columns == 292681L)
stopifnot(clf_features == 292681L)
stopifnot(pgf_probesets == 36353L)
```

The validated outputs carried into the package build are:

- platform name: `miRNA-4_1`;
- geometry: 541 × 541;
- physical features: 292,681;
- PGF probesets: 36,353.

---

## 6. Build the platform package

The PGF and CLF have now passed structural validation. This section combines
them with standard R package metadata and generates the installable package.

### 6.1. Define package metadata

Use valid maintainer information for a distributable build:

```r
package_metadata <- list(
  version = "0.1.1",
  license = "Artistic-2.0",
  author = "Your Name",
  email = "your.name@example.org",
  biocViews = "AnnotationData",
  chipName = "miRNA-4_1",
  manufacturer = "Affymetrix",
  url = "https://www.thermofisher.com/",
  genomebuild = "not applicable",
  organism = "multiple species",
  species = "multiple species"
)
```

The array contains probes for multiple species, so the platform package
should not claim a single mouse genome build. Species-specific filtering is
performed later using the annotation CSV.

### 6.2. Create the miRNA platform seed

Combine the validated files from Section 5 with the metadata from Section 6.1:

```r
suppressPackageStartupMessages(library(pdInfoBuilder))

seed <- methods::new(
  "AffyMiRNAPDInfoPkgSeed",
  pgfFile = normalizePath(pgf_file),
  clfFile = normalizePath(clf_file),
  version = package_metadata$version,
  license = package_metadata$license,
  author = package_metadata$author,
  email = package_metadata$email,
  biocViews = package_metadata$biocViews,
  chipName = package_metadata$chipName,
  manufacturer = package_metadata$manufacturer,
  url = package_metadata$url,
  genomebuild = package_metadata$genomebuild,
  organism = package_metadata$organism,
  species = package_metadata$species
)

seed
stopifnot(methods::is(seed, "AffyMiRNAPDInfoPkgSeed"))
```

The seed is the complete build specification used in the next step.

### 6.3. Generate the package source

Create a local, Git-ignored build directory and pass the seed to
`pdInfoBuilder`:

```r
build_root <- file.path("build", "pd_mirna_4_1")
dir.create(build_root, recursive = TRUE, showWarnings = FALSE)

pdInfoBuilder::makePdInfoPackage(
  seed,
  destDir = build_root,
  quiet = FALSE
)
```

The command should generate:

```text
build/pd_mirna_4_1/pd.mirna.4.1/
├── DESCRIPTION
├── R/
├── data/
└── inst/
    └── extdata/
        └── pd.mirna.4.1.sqlite
```

Define the generated source path for the validation and installation steps:

```r
package_source <- file.path(build_root, "pd.mirna.4.1")

stopifnot(dir.exists(package_source))
stopifnot(file.exists(file.path(package_source, "DESCRIPTION")))

list.files(package_source)
read.dcf(file.path(package_source, "DESCRIPTION"))
```

---

## 7. Validate the generated package database

Section 6 generated an R package source and an SQLite database. Validate that
database before installing the package.

```r
generated_sqlite <- file.path(
  package_source,
  "inst",
  "extdata",
  "pd.mirna.4.1.sqlite"
)

stopifnot(file.exists(generated_sqlite))
stopifnot(file.info(generated_sqlite)$size > 0)

connection <- DBI::dbConnect(
  RSQLite::SQLite(),
  generated_sqlite
)

tables <- DBI::dbListTables(connection)

feature_set_count <- DBI::dbGetQuery(
  connection,
  "SELECT COUNT(*) AS n FROM featureSet"
)$n

pm_feature_count <- DBI::dbGetQuery(
  connection,
  "SELECT COUNT(*) AS n FROM pmfeature"
)$n

type_count <- DBI::dbGetQuery(
  connection,
  "SELECT COUNT(*) AS n FROM type_dict"
)$n

DBI::dbDisconnect(connection)

database_summary <- c(
  featureSet = feature_set_count,
  pmfeature = pm_feature_count,
  type_dict = type_count
)

tables
database_summary

stopifnot(all(c("featureSet", "pmfeature", "type_dict") %in% tables))
stopifnot(feature_set_count == pgf_probesets)
stopifnot(feature_set_count == 36353L)
stopifnot(pm_feature_count == 346085L)
stopifnot(type_count == 5L)
```

This step connects the generated package back to the manufacturer design:
the SQLite `featureSet` count must equal the PGF probeset count validated in
Section 5.3.

---

## 8. Install and verify the package

The package source has now passed database validation. Install that exact
source into the active environment established in Section 4:

```r
if (requireNamespace("renv", quietly = TRUE)) {
  renv::install(package_source)
} else {
  install.packages(
    package_source,
    repos = NULL,
    type = "source"
  )
}
```

Alternatively:

```bash
R CMD INSTALL build/pd_mirna_4_1/pd.mirna.4.1
```

Confirm that the active R session finds both the package and its database:

```r
stopifnot(requireNamespace("pd.mirna.4.1", quietly = TRUE))

packageVersion("pd.mirna.4.1")
find.package("pd.mirna.4.1")

installed_sqlite <- system.file(
  "extdata",
  "pd.mirna.4.1.sqlite",
  package = "pd.mirna.4.1"
)

stopifnot(nzchar(installed_sqlite))
stopifnot(file.exists(installed_sqlite))

installed_sqlite
```

Do not copy the generated source folder directly into an R library. Installing
it properly ensures that R registers the package and that `oligo` can find it.

---

## 9. Validate the package with CEL files

The structural package is now installed. The next test uses real CEL data to
confirm that the package name and physical geometry agree with the experiment.

Define the CEL paths:

```r
cel_dir <- file.path(
  "data",
  "raw",
  "mirna_150001",
  "CEL Files"
)

cel_files <- list.files(
  cel_dir,
  pattern = "[.]CEL$",
  full.names = TRUE,
  ignore.case = TRUE
)

stopifnot(length(cel_files) > 0L)
```

Read one CEL file first:

```r
test_cel <- oligo::read.celfiles(cel_files[[1]])
test_values <- Biobase::exprs(test_cel)

test_cel
dim(test_values)
table(is.finite(test_values))

stopifnot(nrow(test_values) == clf_features)
stopifnot(nrow(test_values) == 292681L)
stopifnot(all(is.finite(test_values)))
stopifnot(Biobase::annotation(test_cel) == "pd.mirna.4.1")
```

If that succeeds, load the complete experiment:

```r
raw_cel <- oligo::read.celfiles(cel_files)
raw_values <- Biobase::exprs(raw_cel)

raw_cel
dim(raw_values)
table(is.finite(raw_values))

stopifnot(nrow(raw_values) == clf_features)
stopifnot(ncol(raw_values) == length(cel_files))
stopifnot(all(is.finite(raw_values)))
stopifnot(Biobase::annotation(raw_cel) == "pd.mirna.4.1")
```

Successful CEL loading confirms that the installed package can map the real
array features. It does not guarantee that every normalization or
background-correction calculation will be numerically stable.

---

## 10. Import the biological annotation

At this point, the platform package has supplied stable probeset identifiers.
The NetAffx CSV introduced in Section 2.3 can now add biological information.

The manufacturer CSV begins with comment lines. Import it with
`comment.char = "#"`:

```r
stopifnot(file.exists(annotation_file))

annotation <- read.csv(
  annotation_file,
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot("Probe Set ID" %in% names(annotation))

annotation$`Probe Set ID` <- as.character(
  annotation$`Probe Set ID`
)

dim(annotation)
names(annotation)
head(
  annotation[
    ,
    c("Probe Set ID", "Probe Set Name", "Sequence Type")
  ]
)

stopifnot(nrow(annotation) == pgf_probesets)
stopifnot(!anyDuplicated(annotation$`Probe Set ID`))
```

After creating an expression table, join by probeset ID:

```r
expression_table$ProbeSetName <- as.character(
  expression_table$ProbeSetName
)

annotated_expression <- merge(
  expression_table,
  annotation,
  by.x = "ProbeSetName",
  by.y = "Probe Set ID",
  all.x = TRUE,
  sort = FALSE
)
```

The platform package and annotation CSV have different jobs:

- `pd.mirna.4.1` maps CEL features to probesets;
- the NetAffx CSV maps probesets to biological labels.

---

## 11. Expected outputs and checkpoints

| Stage | Output or checkpoint | Expected value |
|---|---|---:|
| Environment | Required build packages available | All `TRUE` |
| Design validation | PGF/CLF platform | `miRNA-4_1` |
| Design validation | Array geometry | 541 × 541 |
| Design validation | Physical features | 292,681 |
| Design validation | PGF probesets | 36,353 |
| Package build | Generated source | `build/pd_mirna_4_1/pd.mirna.4.1/` |
| Database validation | SQLite `featureSet` | 36,353 |
| Database validation | SQLite `pmfeature` | 346,085 |
| Database validation | SQLite `type_dict` | 5 |
| Installation | Package name | `pd.mirna.4.1` |
| Installation | Local version used here | 0.1.1 |
| CEL validation | Features per array | 292,681 |
| Annotation | Unique `Probe Set ID` rows | 36,353 |

---

## 12. Troubleshooting

Follow the same order as the build. Preserve error messages,
`sessionInfo()`, and diagnostic output before changing preprocessing choices.

### 12.1. External-file problems

#### 12.1.1. The miRNA 4.1 files are not visible

Search Thermo Fisher’s library-file page for **MiRNA 4.1 Analysis Files for
Expression Console and AGCC**. If unavailable, request that legacy analysis
bundle from Thermo Fisher support. Do not substitute only the TAC
configuration file.

#### 12.1.2. PGF and CLF checks disagree

Stop before package construction. Both files must identify `miRNA-4_1` and
the 541 × 541 design. Re-download a matching pair from the same analysis
bundle.

### 12.2. Package-build problems

#### 12.2.1. An unexpected package name is generated

Confirm:

```r
package_metadata$chipName
```

It must be `miRNA-4_1`, which `pdInfoBuilder` converts to
`pd.mirna.4.1`.

#### 12.2.2. The build directory contains an earlier package

Build in a new empty directory. Do not overwrite an installed package or
build inside `resources/affymetrix_library_files/`.

#### 12.2.3. SQLite counts do not match

Do not install the package. Confirm the PGF/CLF pair, create a clean build
directory, and rerun `makePdInfoPackage()`.

### 12.3. R environment and installation problems

#### 12.3.1. R cannot find `pd.mirna.4.1`

```r
R.version.string
.libPaths()
find.package("pd.mirna.4.1", quiet = TRUE)

if (requireNamespace("renv", quietly = TRUE)) {
  renv::status()
}
```

Install and run the package with the same R environment. Do not permanently
add another R version’s package library to `.libPaths()`.

#### 12.3.2. `read.celfiles()` attempts to download the package

The package is absent from the active library or failed to load:

```r
requireNamespace("pd.mirna.4.1", quietly = TRUE)
system.file(
  "extdata",
  "pd.mirna.4.1.sqlite",
  package = "pd.mirna.4.1"
)
```

Return to Section 8 and reinstall the validated source into the active
environment.

### 12.4. CEL-import problems

#### 12.4.1. CEL files identify another platform

```r
affxparser::readCelHeader(cel_files[[1]])
```

Do not force miRNA 4.1 design files onto CEL files from another platform.

#### 12.4.2. Raw intensities are non-finite

```r
raw_values <- Biobase::exprs(raw_cel)

c(
  NA = sum(is.na(raw_values)),
  NaN = sum(is.nan(raw_values)),
  Inf = sum(is.infinite(raw_values))
)
```

Investigate file integrity, transfer problems, and array-level QC before
normalization.

### 12.5. RMA-normalization problems

#### 12.5.1. RMA fails or produces non-finite values

Confirm the result:

```r
rma_result <- oligo::rma(raw_cel)
rma_expression <- Biobase::exprs(rma_result)

c(
  NA = sum(is.na(rma_expression)),
  NaN = sum(is.nan(rma_expression)),
  Inf = sum(is.infinite(rma_expression))
)
```

Possible causes include:

- incompatible or incorrectly generated platform mappings;
- damaged or unusual CEL data;
- numerical failure during background correction;
- one or more arrays whose intensity distributions cause the combined RMA
  calculation to fail.

If the design and raw-value checks pass, inspect array-level distributions and
run controlled diagnostic subsets to determine whether the failure follows
specific arrays. Keep exclusions diagnostic until metadata and QC support
them.

Evidence-dependent options include:

1. rebuild and revalidate the platform package;
2. obtain replacement CEL files if integrity is questionable;
3. exclude clearly identified problematic arrays with a documented rationale;
4. use manufacturer-processed expression values;
5. test RMA without background correction as a labeled sensitivity analysis.

Keep alternative expression matrices separate and record the preprocessing
source in every downstream analysis.

---

## 13. References

- [Thermo Fisher GeneChip Array Library Files](https://www.thermofisher.com/br/en/home/life-science/microarray-analysis/microarray-data-analysis/genechip-array-library-files.html)
- [Bioconductor `pdInfoBuilder` package](https://bioconductor.org/packages/release/bioc/html/pdInfoBuilder.html)
- [Bioconductor `pdInfoBuilder` manual](https://bioconductor.org/packages/release/bioc/manuals/pdInfoBuilder/man/pdInfoBuilder.pdf)
- [Bioconductor `oligo` package](https://bioconductor.org/packages/release/bioc/html/oligo.html)
- [Bioconductor `oligo` user guide](https://bioconductor.org/packages/release/bioc/vignettes/oligo/inst/doc/oug.pdf)
- [Bioconductor support: Affymetrix miRNA 4.1 and `pd.mirna.4.1`](https://support.bioconductor.org/p/96882/)
