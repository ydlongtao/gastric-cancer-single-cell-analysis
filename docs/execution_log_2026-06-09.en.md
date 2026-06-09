# Execution Log

Date: 2026-06-09

## 1. Project Setup

- Workspace: `/Users/huangfulongtao/Desktop/胃癌单细胞项目`
- Project name for public use: `Gastric Cancer Single-Cell Analysis`
- Created `README.md` and `AGENTS.md`.
- Project scope: reproducible gastric cancer single-cell transcriptomics, spatial transcriptomics, and clinical metadata analysis.

## 2. Environment Checks

- Default Python: `/opt/homebrew/bin/python3`, Python 3.14.5, without Scanpy.
- R: `/usr/local/bin/R`, R 4.6.0.
- Installed R packages: `Seurat` 5.5.0, `SeuratObject` 5.4.0, `scDesign3` 1.10.0.
- micromamba environment: `biomni_e1`.
- Scanpy environment: `/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python`.
- Installed Python packages: `scanpy` 1.11.5, `anndata` 0.12.16, `numpy` 2.1.0, `pandas` 2.3.3.

## 3. Rust Setup

- Installed Rust through Homebrew `rustup`.
- Added `/opt/homebrew/opt/rustup/bin` to `~/.zshrc`.
- Verified `rustup` 1.29.0, `rustc` 1.96.0, and `cargo` 1.96.0.
- Created and tested `scripts/rust/rust_check/`.

## 4. GSE234129 Download and Conversion

- Created `Datasets/`, `Datasets/GSE234129/`, and `Datasets/clinical_records/`.
- Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234129
- GEO supplementary FTP: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234129/suppl/`
- Downloaded locally:
  - `GSE234129_barcodes.tsv.gz`
  - `GSE234129_count_matrix.mtx.gz`
  - `GSE234129_features.tsv.gz`
  - `GSE234129_meta.tsv.gz`
- Created SHA256 checksum records.
- Converted to:
  - `Datasets/GSE234129/processed/GSE234129_scanpy.h5ad`
  - `Datasets/GSE234129/processed/GSE234129_seurat.rds`
- Raw and processed matrix files are local-only and ignored by Git.

## 5. Scanpy Analysis

- User chose Scanpy as the main workflow and global clustering as the first subpopulation analysis scope.
- Added scripts:
  - `scripts/python/01_gse234129_qc.py`
  - `scripts/python/02_gse234129_clustering.py`
  - `scripts/python/03_gse234129_marker_annotation.py`
- Installed missing `leidenalg` into the `biomni_e1` environment.

QC:

- Before QC: 19,488 cells.
- After QC: 19,144 cells.
- Retained all 6 patients and all 17 samples.
- Detected 13 mitochondrial genes.

Clustering:

- 3,000 highly variable genes.
- PCA: 50 components.
- Neighbors: `n_neighbors=15`, `n_pcs=30`.
- Leiden clusters:
  - `0.3`: 14 clusters
  - `0.5`: 20 clusters
  - `0.8`: 27 clusters
  - `1.0`: 28 clusters

Marker and annotation:

- Used Wilcoxon marker ranking on `leiden_0.5`.
- Generated cluster-celltype crosstab and classic marker scores.
- Added `global_annotation` to the annotated object.

## 6. scop Plotting and Cairo PDF Support

- Downloaded `mengxu98/scop` to `external/scop/` using a GitHub zip archive.
- Installed dependencies using `pak`.
- Initial `scop` source build failed because R expected a missing `/opt/gfortran` runtime.
- Installed Homebrew `gcc` and patched local `external/scop/src/Makevars` to link against Homebrew Fortran libraries.
- Installed `scop` 0.8.9 successfully.
- Added:
  - `scripts/python/export_gse234129_scop_inputs.py`
  - `scripts/R/04_gse234129_scop_plots.R`
- Used `CellDimPlot`, `FeatureDimPlot`, `FeatureStatPlot`, and `CellStatPlot`.
- Generated QC, UMAP, mitochondrial percentage, and composition figures.
- Updated `AGENTS.md` so final figures should be exported as paired PDF and TIF files.
- Installed XQuartz 2.8.5 to satisfy R Cairo PDF dependencies.
- Verified `grDevices::cairo_pdf()` works and regenerated PDF/TIF outputs.

## 7. GitHub Publication Preparation

- Target public repository name: `gastric-cancer-single-cell-analysis`.
- Confirmed GitHub CLI login as `ydlongtao`.
- Updated `.gitignore` to exclude raw matrices, processed objects, local result objects, third-party source caches, `.DS_Store`, `.Rhistory`, and `Rplots.pdf`.
- Public repository should contain scripts, documentation, logs, download links, compact summary tables, and final PDF/TIF figures.
