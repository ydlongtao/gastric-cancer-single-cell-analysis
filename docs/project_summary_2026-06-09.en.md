# Gastric Cancer Single-Cell Analysis Summary

Updated: 2026-06-09

## Project

This project provides a reproducible analysis workflow for gastric cancer single-cell transcriptomics, with planned support for spatial transcriptomics and de-identified clinical metadata.

The current first dataset is public GEO dataset GSE234129. The workflow has completed data download, object conversion, Scanpy-based QC/clustering/annotation, scop-based plotting, and GitHub publication preparation.

## Environment

- Python/Scanpy: `/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python`
- Scanpy: 1.11.5
- AnnData: 0.12.16
- Seurat: 5.5.0
- scDesign3: 1.10.0
- scop: 0.8.9
- Rust: `rustc` 1.96.0 and `cargo` 1.96.0
- Cairo PDF support: enabled after installing XQuartz 2.8.5

## Data

GSE234129:

- GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234129
- Supplementary FTP: `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234129/suppl/`
- Scale: 19,488 cells and 27,176 genes/features.
- Metadata: `patient`, `sample`, `celltype`.
- Public GitHub policy: keep raw and processed matrices local only; publish links, checksums, and conversion scripts.

## Implemented Workflows

Scanpy:

- `scripts/python/01_gse234129_qc.py`
- `scripts/python/02_gse234129_clustering.py`
- `scripts/python/03_gse234129_marker_annotation.py`

scop plotting:

- `scripts/python/export_gse234129_scop_inputs.py`
- `scripts/R/04_gse234129_scop_plots.R`

Conversion scripts:

- `Datasets/GSE234129/scripts/make_scanpy_h5ad.py`
- `Datasets/GSE234129/scripts/make_seurat_rds.R`

## Main Results

QC:

- Before QC: 19,488 cells.
- After QC: 19,144 cells.
- Retained all 6 patients and 17 samples.

Clustering:

- Main clustering key: `leiden_0.5`.
- Main cluster count: 20 clusters.
- Additional resolutions: 14, 27, and 28 clusters at `0.3`, `0.8`, and `1.0`.

Annotation:

- Wilcoxon marker ranking for all `leiden_0.5` clusters.
- Cross-checked clusters against existing `celltype` labels and classic marker sets.

Figures:

- QC violin by sample.
- UMAP by Leiden cluster.
- UMAP by global annotation.
- UMAP by sample.
- UMAP by mitochondrial percentage.
- Annotation composition by sample.
- Final figure format: paired PDF and 600 dpi TIF.

## Next Suggested Tasks

- Review UMAPs for patient/sample-driven batch effects.
- Manually refine cluster labels from `GSE234129_cluster_annotation.tsv`.
- Perform second-level clustering for immune, myeloid, B/plasma, or epithelial/tumor-related compartments.
- Add de-identified clinical metadata and test composition or marker differences across clinical groups.
