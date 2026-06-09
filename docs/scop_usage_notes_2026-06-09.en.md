# scop Installation and Usage Notes

Date: 2026-06-09

## Package

- GitHub: https://github.com/mengxu98/scop
- Local source cache: `external/scop/`
- Installed version: `scop` 0.8.9
- R requirement: R >= 4.1.0

## Installation Notes

1. `git clone` was slow, so the GitHub zip archive was downloaded and unpacked to `external/scop/`.
2. `pak::local_install("external/scop", dependencies = TRUE)` installed most dependencies.
3. The first `scop` source build failed because the local R installation expected Fortran runtime libraries under `/opt/gfortran`.
4. Homebrew `gcc` was installed to provide `gfortran` and Fortran runtime libraries.
5. Local `external/scop/src/Makevars` was patched to link against Homebrew gcc libraries.
6. `R CMD INSTALL external/scop` then installed `scop` successfully.
7. XQuartz 2.8.5 was installed to provide `/opt/X11` libraries required by R Cairo PDF output.
8. `grDevices::cairo_pdf()` was verified and is now used by the scop plotting script.

Verification:

```bash
Rscript -e 'library(scop); packageVersion("scop")'
```

## Useful scop Functions Learned

- `CellDimPlot()`: low-dimensional cell plots such as UMAP.
- `FeatureDimPlot()`: continuous feature or metadata visualization on UMAP.
- `FeatureStatPlot()`: violin/box/bar/dot plots for features or metadata.
- `CellStatPlot()`: cell count or composition plots from metadata.
- `RunCellQC()`: built-in scop QC workflow, not used for this first pass because Scanpy QC was already completed.

## Project Inputs

Scanpy annotated object:

```text
results/GSE234129/objects/GSE234129_annotated.h5ad
```

Exported metadata for R/scop:

```text
results/GSE234129/tables/GSE234129_scop_plot_metadata.tsv
```

Temporary Seurat plotting object:

```text
results/GSE234129/objects/GSE234129_scop_plotting_seurat.rds
```

## Reproduction

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/export_gse234129_scop_inputs.py
Rscript scripts/R/04_gse234129_scop_plots.R
```

## Figure Outputs

Directory:

```text
results/GSE234129/figures/scop/
```

Each basename is exported as `.pdf` and `.tif`:

- `GSE234129_scop_qc_violin_by_sample`
- `GSE234129_scop_umap_leiden05`
- `GSE234129_scop_umap_global_annotation`
- `GSE234129_scop_umap_sample`
- `GSE234129_scop_umap_pct_counts_mt`
- `GSE234129_scop_annotation_composition_by_sample`

PDF is the editable vector format. TIF is exported at 600 dpi for publication-style bitmap output. PNG files are preview-only.
