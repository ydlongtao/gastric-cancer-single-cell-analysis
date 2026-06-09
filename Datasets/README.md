# Datasets

本目录用于存放项目相关数据，包括单细胞转录组数据、空间转录组数据和病历/临床资料。

## 目录约定

- `GSE234129/`：GEO 数据集 GSE234129 的原始文件、转换脚本和可分析对象。
- `clinical_records/`：病历资料、临床分组、样本对应表和数据字典。

## 当前数据集

### GSE234129

来源：

```text
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234129
```

已整理文件：

- Seurat：`GSE234129/processed/GSE234129_seurat.rds`
- Scanpy/AnnData：`GSE234129/processed/GSE234129_scanpy.h5ad`

数据规模：

- 19,488 cells
- 27,176 genes/features
- 6 patients
- 17 samples
- 62 annotated cell types

原始文件校验值记录在：

```text
GSE234129/raw/SHA256SUMS.txt
```
