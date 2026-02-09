# MM Myeloid High Antigen Presentation & High Immunosuppressive State Signature (scRNA → bulk score → survival validation)

This repository derives a myeloid functional-state gene signature of the High Antigen Presentation (AP) & High Immunosuppressive (SUP) State Signature from multiple myeloma (MM) scRNA-seq, builds a bulk RNA-seq scoring model (GSVA) with an outcome-driven cutoff, and validates the model in an independent cohort.

## Overview

Pipeline steps:

1. Download GEO datasets (scRNA discovery/validation and bulk validation).
2. Convert GEO raw files into Seurat objects.
3. In discovery scRNA (GSE124310), integrate samples (Harmony), subset myeloid cells, define AP/SUP functional states, and identify AP_high_Sup_high marker genes.
4. In bulk cohort (GSE136337), compute GSVA state score, learn a risk cutoff (maxstat), perform survival analysis, and save a reusable risk model.
5. Validate the saved risk model in an independent cohort (TCGA-LAML) using the fixed cutoff.

All figures are saved to `../results/figures/` as:
- SVG (vector)
- PNG (600 dpi)

Figure styling and saving are controlled in `utils.R`.

## Repository layout

After running the R scripts in the `scripts/` directory, you should expected the below directory structure relative to the `scripts/` directory (or where you run the scripts):

- `../data/geo/`  
  GEO downloads (supplementary files and series matrices)

- `../data/untared_geo/`  
  Untar/unzipped GEO raw files

- `../data/seurat/`  
  Saved Seurat object lists (`*.rds`)

- `../data/TCGA_LAML/`  
  TCGA download cache (created by TCGAbiolinks)

- `../results/`  
  Output tables and models (signature marker lists, risk model)

- `../results/figures/`  
  All figures (SVG + 600 dpi PNG)

## Requirements

- R (recommended: R >= 4.3)
- `renv` for dependency management
- Internet access for GEO downloads and TCGA downloads (TCGAbiolinks)

Key packages (installed via `renv.lock`):
- Seurat, harmony, SingleR, celldex
- GEOquery
- GSVA, survival, survminer
- TCGAbiolinks, SummarizedExperiment, edgeR
- ggplot2, ggrepel, ggpubr, svglite

## Setup (recommended)

From the project root:

```r
install.packages("renv")
renv::restore()
```

This installs the versions pinned in `renv.lock`.

## Scripts

### `01_download_geo.R`

Downloads:
- scRNA discovery: `GSE124310` supplementary files
- scRNA validation: `GSE163278` supplementary files
- bulk validation: `GSE136337` series matrix (saved as RDS)

Outputs:
- `../data/geo/GSE124310/*`
- `../data/geo/GSE163278/*`
- `../data/geo/GSE136337/GSE136337_GSEMatrix.rds`

Run:
```r
source("01_download_geo.R")
```

### `02_scRNA_geo2seurat.R`

Converts GEO raw files into Seurat object lists:

- For `GSE124310`: untar and read 10X matrices, save `GSE124310.seurat.list.rds`
- For `GSE163278`: gunzip/rename matrices to fit `Read10X`, save `GSE163278.seurat.list.rds`

Outputs:
- `../data/seurat/GSE124310.seurat.list.rds`
- `../data/seurat/GSE163278.seurat.list.rds`

Run:
```r
source("02_scRNA_geo2seurat.R")
```

### `03_mm_state_signature.R`

Discovery in scRNA (GSE124310):

- QC filtering (features, mito)
- Harmony integration
- Clustering and UMAP
- SingleR annotation to identify myeloid populations
- Subset myeloid clusters and recluster
- Define functional states based on two module scores:
  - AP (antigen presentation) score from `AP_GENES_REF`
  - SUP (immunosuppressive) score from `SUP_GENES_REF`
- Compute per-stage and per-cluster functional-state composition
- Differential expression to identify AP_high_Sup_high markers
- Export marker tables for bulk scoring

Key outputs:
- `../results/myeloid_APSupHigh_pos_markers.tsv`
- `../results/myeloid_APSupHigh_novel_markers.tsv`
- figures saved under `../results/figures/`

Run:
```r
source("utils.R")
source("03_mm_state_signature.R")
```

### `04_bulk_signature_score.R`

Bulk scoring and survival analysis (GSE136337):

- Load series matrix expression and phenotype data
- Select refined novel markers from `03` output
- Compute GSVA scores:
  - `State_Score`
  - `Myeloid_Abundance` (control set: CD14, CSF1R, CD68, LYZ)
- Determine a risk cutoff with `surv_cutpoint`
- Kaplan–Meier and Cox models (uni/multi), including PH checks
- Save a reusable model object for validation

Key outputs:
- `../results/Myeloid_State_Risk_Model.rds`
- figures saved under `../results/figures/`

Run:
```r
source("utils.R")
source("04_bulk_signature_score.R")
```

### `05_bulk_score_independent_test.R`

Independent validation in TCGA-LAML:

- Download and prepare TCGA-LAML expression counts (STAR counts)
- Convert to log2-CPM and map Ensembl IDs to gene symbols
- Load model from `../results/Myeloid_State_Risk_Model.rds`
- Compute GSVA state scores and apply the fixed cutoff from the model
- Survival curves and Cox models (uni/multi), including PH checks
- figures saved under `../results/figures/`

Run:
```r
source("utils.R")
source("05_bulk_score_independent_test.R")
```

## Figure generation

All scripts are expected to use `save_figure()` from `utils.R` to save:
- `*.svg` via `svglite`
- `*.png` at 600 dpi

If you add new plots, use:
```r
p <- <ggplot or Seurat plot object>
save_figure(p, "../results/figures", "<figure_name>", width = <in>, height = <in>)
```

## Reproducibility notes

- `renv.lock` pins package versions. It is crucial that you use Seurat v4 instead of v5.
- GEO and TCGA downloads depend on remote availability and may change over time.
- Some cohorts may have missing genes; scripts report the percent of model genes covered in the validation dataset.

## Citation / data sources

- GEO: GSE124310, GSE163278, GSE136337
- TCGA: TCGA-LAML (downloaded via TCGAbiolinks)
