# GSE234129

Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234129

Downloaded on 2026-06-09 from the GEO supplementary FTP directory:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234129/suppl/
```

## Raw Files

- `raw/GSE234129_count_matrix.mtx.gz`
- `raw/GSE234129_features.tsv.gz`
- `raw/GSE234129_barcodes.tsv.gz`
- `raw/GSE234129_meta.tsv.gz`
- `raw/SHA256SUMS.txt`

Matrix dimensions:

- genes/features: 27,176
- cells/barcodes: 19,488
- non-zero count entries: 31,652,388

Metadata columns:

- `cell_barcodes`
- `patient`
- `sample`
- `celltype`

## Processed Files

Generate Scanpy/AnnData:

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python Datasets/GSE234129/scripts/make_scanpy_h5ad.py
```

Output:

```text
processed/GSE234129_scanpy.h5ad
```

Generate Seurat:

```bash
Rscript Datasets/GSE234129/scripts/make_seurat_rds.R
```

Output:

```text
processed/GSE234129_seurat.rds
```
