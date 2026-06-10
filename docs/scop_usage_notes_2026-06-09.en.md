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

## GSE234129 Extended scop Analysis

New script:

```bash
Rscript scripts/R/08_gse234129_scop_extended_analysis.R
```

Analyses:

- `GroupHeatmap()` for classic markers, `global_annotation` top markers, and Leiden 0.5 top markers.
- `FeatureHeatmap()` for sampled classic marker expression by `global_annotation`.
- `CellCorHeatmap()` for QC-object comparisons across `leiden_0.5`, `celltype`, `global_annotation`, and `sample`, plus raw-vs-QC comparisons.
- `RunSlingshot()` pseudotime inference for eligible T/NK, Myeloid/TAM, and B/Plasma lineages.
- `PseudotimeProjectionPlot()` and `DynamicHeatmap()` for pseudotime direction and lineage-associated expression dynamics.

Outputs:

```text
results/GSE234129/figures/scop_extended/
results/GSE234129/tables/scop_extended/
results/GSE234129/reports/GSE234129_scop_extended_analysis_report.md
```

Notes:

- The script fixes `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, and `BLIS_NUM_THREADS` to 1 to avoid runaway parallelism.
- `monocle3`, `tradeSeq`, and `palantir` are not currently installed; the formal trajectory workflow uses the installed `slingshot` package.
- When `grouping.var` is used, `GroupHeatmap()` follows `scop` behavior and restricts expression comparisons to `log2fc`.
- B/Plasma Slingshot trajectory outputs are generated successfully; the B/Plasma dynamic heatmap can be locally skipped if several genes fail model fitting because of infinite or missing values.
