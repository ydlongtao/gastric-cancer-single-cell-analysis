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

## 15. scDesign3 四倍模拟执行

日期：2026-06-10

- 新增 scDesign3 输入导出脚本：
  - `scripts/python/export_gse234129_scdesign3_inputs.py`
- 新增 scDesign3 主脚本：
  - `scripts/R/05_gse234129_scdesign3_simulation.R`

执行策略：

- 使用 QC 后 annotated 对象作为正式输入：
  - `results/GSE234129/objects/GSE234129_annotated.h5ad`
- 使用 Seurat counts 对象提供原始计数：
  - `Datasets/GSE234129/processed/GSE234129_seurat.rds`
- 先尝试 `3078` 个 HVG+marker 特征和 `global_annotation + sample`，边际模型拟合过慢。
- 调整为 `200` 个 HVG/marker 特征，并将均值公式降为 `global_annotation`。
- 使用 `n_cores=7`，同时在脚本内将 `OMP_NUM_THREADS`、`OPENBLAS_NUM_THREADS`、`MKL_NUM_THREADS`、`VECLIB_MAXIMUM_THREADS` 等底层线程限制为 `1`，避免递归并行占满线程。

最终命令：

```bash
Rscript scripts/R/05_gse234129_scdesign3_simulation.R \
  --refresh-inputs \
  --max-features=200 \
  --n-cores=7 \
  --mu-formula=global_annotation
```

最终日志：

- `results/GSE234129/scdesign3/logs/GSE234129_scdesign3_4x_200features_global_annotation_20260610_092510.log`

耗时：

- 全流程：25 分 48 秒
- 起止时间：2026-06-10 09:25:17 到 09:51:05
- scDesign3 正式模拟段：2026-06-10 09:25:32 到 09:51:05，约 25 分 33 秒

输出：

- `results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_counts.rds`
- `results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_sce.rds`
- `results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_seurat.rds`
- `results/GSE234129/scdesign3/tables/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_metadata.tsv`
- `results/GSE234129/scdesign3/GSE234129_scdesign3_simulation_report.md`
- QC、composition、UMAP、marker dotplot 的 PDF/TIF/PNG 图件。

验证：

- 模拟 counts：`200 × 76576`
- 模拟 SCE：`200 × 76576`
- 模拟 Seurat：`200 × 76576`
- metadata 行数：76,576
- unique barcodes：76,576
- feature set：200
- `cell_count_ok=TRUE`
- `metadata_ok=TRUE`
- `feature_ok=TRUE`

## 16. 原始数据集 scop 绘图

日期：2026-06-10

- 新增原始数据集 scop metadata 导出脚本：
  - `scripts/python/export_gse234129_raw_scop_inputs.py`
- 新增原始数据集 scop 绘图脚本：
  - `scripts/R/06_gse234129_raw_scop_plots.R`

执行命令：

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/export_gse234129_raw_scop_inputs.py
Rscript scripts/R/06_gse234129_raw_scop_plots.R
```

输入：

- `Datasets/GSE234129/processed/GSE234129_scanpy.h5ad`
- `Datasets/GSE234129/processed/GSE234129_seurat.rds`

方法：

- 保留原始未过滤全集 19,488 个细胞。
- 重新计算 QC 指标：`n_genes_by_counts`、`total_counts`、`pct_counts_mt`。
- 在原始全集上计算 HVG、PCA、neighbors、UMAP 和 Leiden 0.5。
- 使用 scop 输出正式 PDF/TIF 图件，PNG 仅作预览。

验证：

- exported metadata：19,488 行，10 列
- Seurat plotting object：`27176 × 19488`
- UMAP cells：19,488
- samples：17
- celltypes：62
- Leiden 0.5 clusters：19

输出：

- `results/GSE234129/tables/GSE234129_raw_scop_plot_metadata.tsv`
- `results/GSE234129/objects/GSE234129_raw_scop_plotting_seurat.rds`
- `results/GSE234129/figures/scop_raw/GSE234129_raw_scop_umap_leiden05.pdf`
- `results/GSE234129/figures/scop_raw/GSE234129_raw_scop_umap_celltype.pdf`
- `results/GSE234129/figures/scop_raw/GSE234129_raw_scop_umap_sample.pdf`
- `results/GSE234129/figures/scop_raw/GSE234129_raw_scop_qc_violin_by_sample.pdf`
- `results/GSE234129/figures/scop_raw/GSE234129_raw_scop_umap_pct_counts_mt.pdf`
- `results/GSE234129/figures/scop_raw/GSE234129_raw_scop_celltype_composition_by_sample.pdf`

每张图同步输出 `.tif` 正式位图版本。
