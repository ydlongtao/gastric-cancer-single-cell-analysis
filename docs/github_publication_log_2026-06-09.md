# GitHub Publication Log

Date: 2026-06-09

## Request

Compress the latest task history into logs, save English versions of README/AGENTS/task logs, rename the project in English, and publish the project to GitHub while avoiding large raw matrix uploads.

## Decisions

- Public project/repository name: `gastric-cancer-single-cell-analysis`.
- Keep raw and processed matrix data local.
- Publish only scripts, documentation, logs, and public download links.
- Use `README.en.md`, `AGENTS.en.md`, and English docs under `docs/` as English copies.

## Large File Policy

Excluded from Git:

- `Datasets/*/raw/`
- `Datasets/*/processed/`
- large `results/` objects and full marker tables
- `external/`
- `.DS_Store`
- `.Rhistory`
- `Rplots.pdf`

Included from `results/`:

- final scop PDF/TIF figures
- compact analysis report
- compact annotation and top-marker summary tables

Raw GSE234129 data can be redownloaded from:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234129/suppl/
```

## GitHub

Target account: `ydlongtao`.

Repository visibility: public.
