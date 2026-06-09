# Gastric Cancer Single-Cell Analysis

This project provides a reproducible workflow for gastric cancer single-cell transcriptomics analysis, with room for spatial transcriptomics and clinical metadata integration. The current implementation uses public GEO dataset GSE234129 as a first test dataset and combines Python/Scanpy, R/Seurat/scop, and optional Rust utilities.

## Project Goals

- Perform quality control, normalization, dimensionality reduction, clustering, and cell population annotation for single-cell RNA-seq data.
- Support future spatial transcriptomics analysis and single-cell/spatial integration.
- Keep raw and processed matrix data local, while preserving public download links and reproducible scripts in Git.
- Produce publication-friendly figures in PDF and TIF formats.

## Environment

Checked on 2026-06-09.

Python/Scanpy:

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python
```

Installed packages in that environment:

- `scanpy` 1.11.5
- `anndata` 0.12.16
- `numpy` 2.1.0
- `pandas` 2.3.3

R:

```bash
/usr/local/bin/R
```

Installed R packages:

- `Seurat` 5.5.0
- `SeuratObject` 5.4.0
- `scDesign3` 1.10.0
- `scop` 0.8.9

Rust:

- `rustup` 1.29.0
- `rustc` 1.96.0
- `cargo` 1.96.0

## Dataset

GSE234129 source:

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234129
```

GEO supplementary download directory:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234129/suppl/
```

Raw and processed matrix files are intentionally excluded from Git. Recreate them with:

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python Datasets/GSE234129/scripts/make_scanpy_h5ad.py
Rscript Datasets/GSE234129/scripts/make_seurat_rds.R
```

Dataset scale:

- 19,488 cells
- 27,176 genes/features
- 6 patients
- 17 samples
- 62 existing cell type labels

## Scanpy Workflow

Run from the project root:

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/01_gse234129_qc.py
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/02_gse234129_clustering.py
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/03_gse234129_marker_annotation.py
```

Main outputs, generated locally and ignored by Git:

- `results/GSE234129/objects/GSE234129_qc_filtered.h5ad`
- `results/GSE234129/objects/GSE234129_clustered.h5ad`
- `results/GSE234129/objects/GSE234129_annotated.h5ad`
- `results/GSE234129/tables/GSE234129_cluster_annotation.tsv`
- `results/GSE234129/tables/GSE234129_leiden05_markers.tsv`

## scop Plotting Workflow

Run:

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/export_gse234129_scop_inputs.py
Rscript scripts/R/04_gse234129_scop_plots.R
```

The plotting script writes paired PDF and TIF files under:

```text
results/GSE234129/figures/scop/
```

PNG files are preview-only and are not the preferred deliverable format.

## Git/Data Policy

- Do not commit raw matrices, processed `.h5ad`/`.rds` objects, local result objects, or large generated figures.
- Keep download links, checksums, scripts, and analysis logs in Git.
- Store clinical files only after de-identification.
