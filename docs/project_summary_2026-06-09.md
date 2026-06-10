# 胃癌单细胞项目摘要

更新时间：2026-06-09

## 项目定位

本项目用于胃癌单细胞转录组、空间转录组和临床病历资料的整合分析。当前第一阶段以公共单细胞数据集 GSE234129 为起点，采用 Scanpy 主流程完成数据质控、降维聚类和全局细胞亚群鉴定。

## 当前环境

- Python/Scanpy 环境：`/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python`
- Scanpy：`1.11.5`
- AnnData：`0.12.16`
- Seurat：`5.5.0`
- scDesign3：`1.10.0`
- Rust：已安装 stable 工具链
  - `rustup 1.29.0`
  - `rustc 1.96.0`
  - `cargo 1.96.0`

注意：默认 `/opt/homebrew/bin/python3` 未安装 Scanpy，单细胞 Python 分析请优先使用 `biomni_e1`。

## 数据集

GSE234129 已下载并整理：

- 来源：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234129
- 原始文件目录：`Datasets/GSE234129/raw/`
- Scanpy 对象：`Datasets/GSE234129/processed/GSE234129_scanpy.h5ad`
- Seurat 对象：`Datasets/GSE234129/processed/GSE234129_seurat.rds`
- 病历资料占位目录：`Datasets/clinical_records/`

原始矩阵规模：

- 19,488 cells
- 27,176 genes/features
- metadata columns: `patient`, `sample`, `celltype`
- 6 patients
- 17 samples
- 62 existing celltype labels

## 已实现的 Scanpy 主流程

脚本：

- `scripts/python/01_gse234129_qc.py`
- `scripts/python/02_gse234129_clustering.py`
- `scripts/python/03_gse234129_marker_annotation.py`

运行命令：

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/01_gse234129_qc.py
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/02_gse234129_clustering.py
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/03_gse234129_marker_annotation.py
```

主要输出：

- QC 后对象：`results/GSE234129/objects/GSE234129_qc_filtered.h5ad`
- 聚类对象：`results/GSE234129/objects/GSE234129_clustered.h5ad`
- 注释对象：`results/GSE234129/objects/GSE234129_annotated.h5ad`
- 分析报告：`results/GSE234129/GSE234129_analysis_report.md`
- 主注释表：`results/GSE234129/tables/GSE234129_cluster_annotation.tsv`
- marker 表：`results/GSE234129/tables/GSE234129_leiden05_markers.tsv`
- UMAP 图：`results/GSE234129/figures/`
- QC 图和保留率表：`results/GSE234129/qc/`

## 分析结果摘要

QC 阈值：

- `n_genes_by_counts >= 300`
- `n_genes_by_counts <= 6000`
- `total_counts >= 500`
- `pct_counts_mt <= 15`

QC 结果：

- QC 前：19,488 cells
- QC 后：19,144 cells
- patient 保留：6/6
- sample 保留：17/17
- mitochondrial genes：13

聚类结果：

- HVGs：3,000
- PCA：50 components
- neighbors：`n_neighbors=15`, `n_pcs=30`
- 主聚类：`leiden_0.5`
- `leiden_0.3`: 14 clusters
- `leiden_0.5`: 20 clusters
- `leiden_0.8`: 27 clusters
- `leiden_1.0`: 28 clusters

亚群鉴定：

- 对 `leiden_0.5` 的 20 个 cluster 进行 Wilcoxon marker gene ranking。
- 使用现有 `celltype` 标签做交叉表。
- 使用经典 marker set 辅助判定全局注释：
  - T/NK
  - B/Plasma
  - Myeloid/TAM
  - Epithelial/Tumor-related
  - Endothelial/Fibroblast

## 下次任务建议入口

优先查看：

1. `docs/project_summary_2026-06-09.md`
2. `docs/execution_log_2026-06-09.md`
3. `results/GSE234129/GSE234129_analysis_report.md`
4. `results/GSE234129/tables/GSE234129_cluster_annotation.tsv`

下一步可做：

- 检查 UMAP 是否存在 patient/sample 批次驱动。
- 根据 `GSE234129_cluster_annotation.tsv` 人工修订 cluster 注释。
- 对 T/NK、髓系、B/Plasma 或 epithelial/tumor-related 细胞分别做二次重聚类。
- 引入病历资料后做样本组间细胞比例和 marker 差异分析。

## scDesign3 四倍模拟更新

更新时间：2026-06-10

新增脚本：

- `scripts/python/export_gse234129_scdesign3_inputs.py`
- `scripts/R/05_gse234129_scdesign3_simulation.R`

最终成功配置：

- 输入对象：`results/GSE234129/objects/GSE234129_annotated.h5ad`
- 基准细胞数：19,144
- 模拟细胞数：76,576
- 特征基因数：200
- 公式：`global_annotation`
- R worker：7
- 底层 BLAS/OMP 线程：1
- 耗时：25 分 48 秒（2026-06-10 09:25:17 到 09:51:05）

主要输出：

- `results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_counts.rds`
- `results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_sce.rds`
- `results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_seurat.rds`
- `results/GSE234129/scdesign3/GSE234129_scdesign3_simulation_report.md`

验证结果：

- counts、SCE、Seurat 对象均为 `200 × 76576`
- metadata 行数为 76,576
- feature set 行数为 200
- 4 倍模拟目标通过

## 原始数据集 scop 绘图更新

更新时间：2026-06-10

新增脚本：

- `scripts/python/export_gse234129_raw_scop_inputs.py`
- `scripts/R/06_gse234129_raw_scop_plots.R`

流程说明：

- 使用 `Datasets/GSE234129/processed/GSE234129_scanpy.h5ad` 作为原始未过滤输入。
- 不套用 QC 过滤阈值，保留全部 19,488 个细胞。
- 在原始全集上重新计算 QC 指标、HVG、PCA、neighbors、UMAP 和 Leiden 0.5。
- 使用 `Datasets/GSE234129/processed/GSE234129_seurat.rds` 与导出的 UMAP/QC metadata 生成 scop 图件。

验证结果：

- Seurat plotting object：`27176 × 19488`
- metadata 行数：19,488
- UMAP 细胞数：19,488
- sample 数：17
- celltype 标签数：62
- Leiden 0.5 cluster 数：19

图件目录：

- `results/GSE234129/figures/scop_raw/`

已生成图件：

- Leiden 0.5 UMAP
- celltype UMAP
- sample UMAP
- QC violin by sample
- mitochondrial percentage UMAP
- celltype composition by sample
