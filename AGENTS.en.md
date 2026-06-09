# AGENTS.md

This file guides future agents and collaborators working on this project.

## Project Positioning

This is a gastric cancer single-cell and spatial transcriptomics analysis project. The default objective is to produce reproducible scripts, figures, tables, and reports rather than only giving instructions.

## Languages and Environment

- R: Seurat, scDesign3, scop, spatial transcriptomics statistics, and publication-style figures.
- Python: Scanpy, AnnData, format conversion, batch processing, and supplementary analysis.
- Rust: optional high-performance command-line utilities.

Confirmed local environment:

- Python/Scanpy: `/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python`
- R/Rscript: `/usr/local/bin/R`, `/usr/local/bin/Rscript`
- R packages: `Seurat` 5.5.0, `SeuratObject` 5.4.0, `scDesign3` 1.10.0, `scop` 0.8.9
- Rust: `rustc` 1.96.0, `cargo` 1.96.0, `rustup` 1.29.0

The default `/opt/homebrew/bin/python3` does not have Scanpy installed. Use the `biomni_e1` environment for Python single-cell analysis.

## Working Rules

- Inspect the directory structure and existing files before creating or modifying files.
- Treat raw data as read-only.
- Be careful with large sequencing files, images, matrix objects, and intermediate objects; do not copy or move them casually.
- Keep scripts reproducible and avoid hidden notebook state.
- Use traceable file names that include the dataset, method, grouping, or date.
- Preserve previous results or document why they changed unless the user explicitly asks for cleanup.

## Directory Conventions

- `Datasets/`: public dataset metadata, download notes, conversion scripts, and clinical-data placeholders.
- `Datasets/*/raw/`: local raw data, ignored by Git.
- `Datasets/*/processed/`: local converted objects, ignored by Git.
- `scripts/R/`: R scripts.
- `scripts/python/`: Python scripts.
- `scripts/rust/`: Rust utility source code.
- `results/`: generated local outputs, ignored by Git.
- `docs/`: documentation, run logs, summaries, and data dictionaries.

## Environment Checks

Scanpy:

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python - <<'PY'
import scanpy as sc
print(sc.__version__)
PY
```

Seurat/scDesign3/scop:

```bash
Rscript -e 'library(Seurat); library(scDesign3); library(scop); sessionInfo()'
```

Rust:

```bash
cargo --version
rustc --version
cargo test --manifest-path scripts/rust/rust_check/Cargo.toml
```

## Code Style

- R: use clear function boundaries and split long workflows into QC, integration, annotation, statistics, and plotting sections.
- Python: prefer `pathlib`, explicit input/output paths, and `argparse` when a script needs parameters.
- Rust: keep tools small and CLI-oriented, and document input/output formats in README files.
- Plotting scripts should fix random seeds when relevant and record key plotting parameters.

## Figure Output Rules

- Prefer publication- and presentation-friendly `PDF` and `TIF/TIFF` outputs.
- For the same figure, save one vector `.pdf` and one high-resolution `.tif`.
- Use at least 300 dpi for TIF; use 600 dpi for publication or detailed layout work.
- Use PNG only for quick previews or web display, not as the primary deliverable.
- Save figures under `results/.../figures/` with descriptive dataset, method, and grouping names.

## Single-Cell and Spatial Analysis Notes

- Record species, reference genome, sample source, group information, and batch information.
- Do not scatter QC thresholds across multiple scripts; record them near the top of the workflow or in a config file.
- Save marker evidence, reference datasets, or manual correction notes for cell type annotation.
- For spatial transcriptomics, record coordinate systems, section IDs, spot/cell filtering rules, and tissue region definitions.
- Keep per-sample QC reports before cross-sample integration.

## Delivery Checklist

When completing a task, report:

- Changed script or document paths.
- Environment and key package versions.
- Generated result paths.
- Verification commands that succeeded.
- Any remaining limitations or user-confirmation needs.
