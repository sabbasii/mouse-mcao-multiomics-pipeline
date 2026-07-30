# Mouse Whole-Blood miRNA–mRNA Microarray Analysis Pipeline

Reproducible R/Bioconductor workflows for analyzing matched mouse whole-blood
miRNA and mRNA microarray profiles during transient middle cerebral artery
occlusion (MCAO), including quality control, normalization, annotation,
differential expression, pathway analysis, and late multi-omics integration.

## Overview

This project analyzes peripheral whole-blood expression during a transient
middle cerebral artery occlusion (MCAO) experiment. The primary paired cohort
contains 43 mice in three groups:

| Group | Description | n |
|---|---|---:|
| Sham | Sham surgery with recovery before collection | 16 |
| MCAO1hr | 60 minutes of occlusion, collected immediately | 16 |
| MCAO3hr | 60 minutes of occlusion, recanalization, and 3 hours of recovery | 11 |

The molecular layers are analyzed independently before late integration. The
primary expression model is:

```text
expression ~ treatment + sex + age
```

The prespecified contrasts are:

- `MCAO1hr_vs_Sham`
- `MCAO3hr_vs_Sham`
- `MCAO3hr_vs_MCAO1hr`

The integrated analysis combines differential-expression evidence with
miRTarBase and TargetScanMouse target evidence and adjusted paired-animal
miRNA–mRNA associations. These analyses prioritize hypotheses; they do not
establish direct regulation, causality, mediation, or cell abundance.

## Public workflow

The complete step-by-step workflow is documented in:

[notebooks/microarray_analysis_pipeline.qmd](notebooks/microarray_analysis_pipeline.qmd)

View the rendered workflow on [GitHub Pages](https://sabbasii.github.io/mouse-mcao-multiomics-pipeline/microarray_analysis_pipeline.html).

Render it with:

```bash
quarto render notebooks/microarray_analysis_pipeline.qmd
```

Private manuscript-preparation files are kept under `notebooks/` locally but
are excluded by the root `.gitignore`.

## Repository structure

```text
scripts/
├── mirna/                 miRNA preprocessing, QC, annotation, and DE
├── mrna/                  mRNA preprocessing, QC, annotation, and DE
├── multiomics/            target evidence, integration, association, and pathways
├── manuscript_figures/    Reproducible publication-figure scripts
└── neutrophil_analysis/   Focused neutrophil-marker analyses

notebooks/                 Public Quarto workflow and supporting images
metadata/                  Placeholder only; local metadata files are not tracked
resources/                 Downloaded/local annotation and pathway resources
data/                      Local raw and processed array data (not tracked)
results/                   Local generated outputs (not tracked)
```

## Analysis sequence

1. Build and validate miRNA and mRNA sample sheets.
2. Normalize and quality-control each molecular layer separately.
3. Annotate mRNA transcript clusters and prepare the unique gene-mapped matrix.
4. Freeze the paired animal manifest.
5. Run treatment, sex, and age-adjusted differential expression.
6. Prepare and validate miRTarBase 10.0 and TargetScanMouse 8.0 evidence.
7. Match target evidence to exact measured miRNA and mRNA features.
8. Perform late integration and adjusted paired-animal association testing.
9. Run full-ranked mRNA GSEA and formal integrated ORA.
10. Run focused manuscript-figure and neutrophil-marker analyses as needed.

Scripts are numbered within their analysis area. The Quarto workflow documents
the corresponding inputs, outputs, and validation checks.

## Requirements

- R 4.6.0 or a compatible recent R release;
- Quarto for rendering the workflow document;
- a working Bioconductor installation;
- the project-local `renv` environment;
- access to the required local array files and annotation resources.

Restore the recorded R environment from the repository root:

```r
install.packages("renv")
renv::restore()
```

## Data and generated outputs

Raw CEL files, vendor CHP files, normalized matrices, pathway databases,
target-resource downloads, QC files, figures, and analysis tables are kept
locally and are excluded from Git. This prevents large or machine-specific
files from entering the public repository.

The mRNA source data are available through GEO accession
[GSE278554](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278554).
The miRNA accession and final analysis-code repository citation will be added
when available.

Resource download instructions, release information, and checksums should be
recorded with the corresponding preparation scripts rather than committing
the downloaded databases themselves.

## Reproducibility notes

- Sample matching is performed using documented identifiers and exact file
   names before biological comparisons are run.
- The primary paired analysis excludes unpaired samples and the single
   MCAO24hr animal from group-level inference.
- miRNA DABG values are used only for independent detection filtering; they are
   not expression measurements or model covariates.
- Results from bulk whole blood cannot distinguish cell-composition changes
   from cell-intrinsic transcriptional regulation.
- Multiple-testing correction and exploratory thresholds are reported
   separately throughout the workflow.

## Citation

Please cite the associated manuscript and the data source when using this
pipeline. A complete citation will be added when the manuscript and public
code release are finalized.

## License

This project is released under the [MIT License](LICENSE).
