# 执行日志

日期：2026-06-09

## 1. 项目初始化

- 工作目录：`/Users/huangfulongtao/Desktop/胃癌单细胞项目`
- 新建项目文档：
  - `README.md`
  - `AGENTS.md`
- 项目用途：使用 R、Python、Rust 及相关程序包分析单细胞转录组和空间转录组数据。

## 2. 环境检查

- 默认 Python：`/opt/homebrew/bin/python3`
  - Python `3.14.5`
  - 未安装 `scanpy`
- R：`/usr/local/bin/R`
  - R `4.6.0`
  - `Seurat 5.5.0`
  - `SeuratObject 5.4.0`
  - `scDesign3 1.10.0`
- micromamba：
  - 可用环境：`biomni_e1`
  - Python 路径：`/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python`
  - `scanpy 1.11.5`
  - `anndata 0.12.16`
  - `numpy 2.1.0`
  - `pandas 2.3.3`

## 3. Rust 安装

- 使用 Homebrew 安装 `rustup`。
- 安装 stable Rust 工具链。
- 写入 `~/.zshrc`：

```bash
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
```

- 验证版本：
  - `rustup 1.29.0`
  - `rustc 1.96.0`
  - `cargo 1.96.0`
  - `rustfmt 1.9.0`
  - `clippy 0.1.96`
- 新建并验证 Rust 测试工程：
  - `scripts/rust/rust_check/`
  - `cargo clippy` 通过
  - `cargo test` 通过
  - `cargo run` 输出 `Hello, world!`

## 4. GSE234129 下载与整理

- 新建数据目录：
  - `Datasets/`
  - `Datasets/GSE234129/raw/`
  - `Datasets/GSE234129/processed/`
  - `Datasets/GSE234129/scripts/`
  - `Datasets/clinical_records/`

- GEO 页面：
  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234129

- GEO supplementary FTP：
  - `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234129/suppl/`

- 下载文件：
  - `GSE234129_barcodes.tsv.gz`
  - `GSE234129_count_matrix.mtx.gz`
  - `GSE234129_features.tsv.gz`
  - `GSE234129_meta.tsv.gz`

- 生成校验文件：
  - `Datasets/GSE234129/raw/SHA256SUMS.txt`

- 原始矩阵检查：
  - genes/features：27,176
  - cells/barcodes：19,488
  - non-zero entries：31,652,388
  - metadata lines：19,489 including header

## 5. 格式转换

- 新增 Scanpy 转换脚本：
  - `Datasets/GSE234129/scripts/make_scanpy_h5ad.py`
- 输出：
  - `Datasets/GSE234129/processed/GSE234129_scanpy.h5ad`
- 验证：
  - AnnData object with `19488 × 27176`
  - `obs`: `patient`, `sample`, `celltype`

- 新增 Seurat 转换脚本：
  - `Datasets/GSE234129/scripts/make_seurat_rds.R`
- 输出：
  - `Datasets/GSE234129/processed/GSE234129_seurat.rds`
- 验证：
  - Seurat object with `27176 features across 19488 samples`
  - Active assay：RNA counts

## 6. 单细胞分析计划确认

- 用户要求：制定并实现单细胞分析计划，至少完成数据质控和细胞亚群鉴定。
- 已确认选择：
  - 主流程：Scanpy
  - 亚群范围：全局亚群

## 7. Scanpy 分析脚本实现

- 新增 QC 脚本：
  - `scripts/python/01_gse234129_qc.py`
- 新增聚类脚本：
  - `scripts/python/02_gse234129_clustering.py`
- 新增 marker/注释脚本：
  - `scripts/python/03_gse234129_marker_annotation.py`

- 安装缺失依赖：
  - `leidenalg 0.12.0`

## 8. QC 执行结果

- 输入：
  - `Datasets/GSE234129/processed/GSE234129_scanpy.h5ad`
- 输出：
  - `results/GSE234129/objects/GSE234129_qc_filtered.h5ad`
  - `results/GSE234129/qc/GSE234129_qc_report.md`
  - QC metric summary 和 retention tables
  - before/after QC histograms

- 阈值：
  - `n_genes_by_counts >= 300`
  - `n_genes_by_counts <= 6000`
  - `total_counts >= 500`
  - `pct_counts_mt <= 15`

- 结果：
  - before QC：19,488 cells
  - after QC：19,144 cells
  - patient 保留：6
  - sample 保留：17
  - mitochondrial genes：13

## 9. 聚类执行结果

- 输入：
  - `results/GSE234129/objects/GSE234129_qc_filtered.h5ad`
- 输出：
  - `results/GSE234129/objects/GSE234129_clustered.h5ad`
  - `results/GSE234129/tables/GSE234129_leiden_cluster_counts.tsv`
  - UMAP figures

- 参数：
  - `normalize_total(target_sum=1e4)`
  - `log1p`
  - HVGs：3,000
  - PCA：50 components
  - neighbors：`n_neighbors=15`, `n_pcs=30`
  - Leiden resolutions：`0.3`, `0.5`, `0.8`, `1.0`

- 结果：
  - `leiden_0.3`: 14 clusters
  - `leiden_0.5`: 20 clusters
  - `leiden_0.8`: 27 clusters
  - `leiden_1.0`: 28 clusters

## 10. Marker 和亚群注释执行结果

- 输入：
  - `results/GSE234129/objects/GSE234129_clustered.h5ad`
- 输出：
  - `results/GSE234129/objects/GSE234129_annotated.h5ad`
  - `results/GSE234129/tables/GSE234129_leiden05_markers.tsv`
  - `results/GSE234129/tables/GSE234129_leiden05_top20_markers.tsv`
  - `results/GSE234129/tables/GSE234129_leiden05_celltype_crosstab.tsv`
  - `results/GSE234129/tables/GSE234129_classic_marker_scores.tsv`
  - `results/GSE234129/tables/GSE234129_cluster_annotation.tsv`
  - `results/GSE234129/GSE234129_analysis_report.md`

- 方法：
  - 对 `leiden_0.5` 使用 Wilcoxon marker ranking。
  - 使用原始 `celltype` 做 cluster-celltype crosstab。
  - 使用经典 marker set 打分辅助注释。

- 验证：
  - marker 表非空，约 543,521 行含 header。
  - cluster annotation 表包含 20 个 cluster。
  - annotated h5ad 中存在 `global_annotation`。

## 11. 复现命令

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/01_gse234129_qc.py
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/02_gse234129_clustering.py
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/03_gse234129_marker_annotation.py
```

## 12. 当前可继续推进的方向

- 人工检查 `results/GSE234129/figures/` 中的 UMAP 图。
- 基于 `GSE234129_cluster_annotation.tsv` 修订 cluster 命名。
- 对免疫细胞、髓系细胞或 epithelial/tumor-related 细胞做二次重聚类。
- 引入 `Datasets/clinical_records/` 中的临床表后，做组间细胞比例和 marker 差异分析。

## 13. scop 绘图与 Cairo PDF 支持

- 安装并学习 GitHub 包 `mengxu98/scop`。
- `git clone` 较慢，改用 GitHub zip 下载到 `external/scop/`。
- `pak::local_install()` 安装大部分依赖成功，但 `scop` 源码构建最初失败，原因是 R 的 Fortran runtime 指向不存在的 `/opt/gfortran`。
- 通过 Homebrew 安装 `gcc`，获得 `gfortran` 和 Fortran runtime。
- 修改本地 `external/scop/src/Makevars`，将 `$(FLIBS)` 替换为 Homebrew gcc 的实际库路径，之后 `R CMD INSTALL external/scop` 成功。
- 使用 `scop::CellDimPlot()`、`scop::FeatureDimPlot()`、`scop::FeatureStatPlot()`、`scop::CellStatPlot()` 对 GSE234129 结果绘图。
- 新增脚本：
  - `scripts/python/export_gse234129_scop_inputs.py`
  - `scripts/R/04_gse234129_scop_plots.R`
- 输出图目录：
  - `results/GSE234129/figures/scop/`
- 已生成 QC、UMAP、线粒体比例和细胞组成图。
- 根据用户要求，在 `AGENTS.md` 中新增图件输出规则：正式图件优先输出 `.pdf` 和 `.tif`，PNG 仅作预览。
- 首次使用 `cairo_pdf` 时 R 报缺失 `/opt/X11` 依赖库。
- 安装 XQuartz 2.8.5 后，`grDevices::cairo_pdf()` 已验证可用。
- 已重新用 `cairo_pdf` 输出 PDF，并同步输出 600 dpi TIF。

## 14. GitHub 公开同步准备

- 用户要求将项目改为英文名称并公开同步到 GitHub。
- 计划使用英文仓库名：`gastric-cancer-single-cell-analysis`。
- 已确认 GitHub CLI 登录账户：`ydlongtao`。
- 已扩展 `.gitignore`，避免上传：
  - 原始矩阵和压缩矩阵文件
  - `.h5ad`、`.rds` 等本地处理对象
  - `results/` 下的大型分析对象和完整 marker 大表
  - `external/` 下的第三方源码缓存和构建产物
  - `.DS_Store`、`.Rhistory`、`Rplots.pdf` 等本地临时文件
- GitHub 中保留轻量级公开结果：
  - scop 输出的 `.pdf` / `.tif` 正式图件
  - 简短分析报告
  - cluster 注释表和 top marker 摘要表
- 原始矩阵数据在公开仓库中仅保留 GEO/GEO FTP 下载链接和转换脚本。
