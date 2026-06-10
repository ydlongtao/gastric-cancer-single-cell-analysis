# scop 安装与使用笔记

日期：2026-06-09

## 包信息

- GitHub: https://github.com/mengxu98/scop
- 本地源码目录：`external/scop/`
- 安装版本：`scop 0.8.9`
- R 要求：R >= 4.1.0

## 安装过程

1. `git clone` 传输较慢，改用 GitHub zip 包下载并解压到 `external/scop/`。
2. 使用 `pak::local_install("external/scop", dependencies = TRUE)` 安装依赖。
3. 首次构建 `scop` 失败，原因是系统 R 的 Fortran runtime 配置指向不存在的 `/opt/gfortran`。
4. 通过 Homebrew 安装 `gcc`，获得 `gfortran` 和 Fortran runtime。
5. 修改本地源码 `external/scop/src/Makevars`，将 `$(FLIBS)` 替换为 Homebrew gcc 的实际库路径。
6. 使用 `R CMD INSTALL external/scop` 安装成功。
7. 通过 Homebrew/XQuartz pkg 安装 XQuartz 2.8.5，补齐 `/opt/X11` 下的 Cairo PDF 依赖库；`grDevices::cairo_pdf()` 已验证可用。

验证命令：

```bash
Rscript -e 'library(scop); packageVersion("scop")'
```

## 本次学习到的常用函数

- `CellDimPlot()`：基于 Seurat reduction 绘制 UMAP/tSNE/PCA 等低维图。
- `FeatureDimPlot()`：在低维坐标上展示 feature 或 metadata 数值。
- `FeatureStatPlot()`：按分组绘制 feature/metadata 的 violin、box、bar、dot 等统计图。
- `CellStatPlot()`：按细胞 metadata 统计比例或数量，可用于细胞组成图。
- `RunCellQC()`：scop 内置 QC 流程；本项目当前使用已完成的 Scanpy QC 字段作为输入绘图。

## 本项目绘图输入

Scanpy 结果对象：

```text
results/GSE234129/objects/GSE234129_annotated.h5ad
```

导出的 scop 绘图 metadata：

```text
results/GSE234129/tables/GSE234129_scop_plot_metadata.tsv
```

用于 scop 绘图的 Seurat 对象：

```text
results/GSE234129/objects/GSE234129_scop_plotting_seurat.rds
```

## 复现命令

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/export_gse234129_scop_inputs.py
Rscript scripts/R/04_gse234129_scop_plots.R
```

## 输出图片

目录：

```text
results/GSE234129/figures/scop/
```

正式输出格式：

- `.pdf`：矢量图，优先用于汇报、排版和后续编辑。
- `.tif`：600 dpi 高分辨率位图，优先用于投稿或需要位图交付的场景。
- `.png`：仅保留为早期快速预览，不作为主要交付格式。

图片 basename，每个 basename 均输出 `.pdf` 和 `.tif`：

- `GSE234129_scop_qc_violin_by_sample`
- `GSE234129_scop_umap_leiden05`
- `GSE234129_scop_umap_global_annotation`
- `GSE234129_scop_umap_sample`
- `GSE234129_scop_umap_pct_counts_mt`
- `GSE234129_scop_annotation_composition_by_sample`

## 注意事项

- 当前 `FeatureStatPlot()` 绘制 QC metadata 时会提示 `Layer 'data' is empty`，因为输入 Seurat 对象只保留 counts layer；图像经检查正常显示 QC metadata 分布。
- 如果后续要用 `scop::RunCellQC()` 或其他深度流程，建议基于原始 Seurat 对象先运行 `NormalizeData()` 或完整 Seurat/scop 标准化流程。

## GSE234129 扩展 scop 分析

新增脚本：

```bash
Rscript scripts/R/08_gse234129_scop_extended_analysis.R
```

分析内容：

- `GroupHeatmap()`：经典 marker、`global_annotation` top markers、Leiden 0.5 top markers。
- `FeatureHeatmap()`：按 `global_annotation` 抽样展示经典 marker 表达。
- `CellCorHeatmap()`：比较 QC 对象内部的 `leiden_0.5`、`celltype`、`global_annotation`、`sample`，并补充 raw-vs-QC 对照。
- `RunSlingshot()`：对 T/NK、Myeloid/TAM、B/Plasma 三个满足规模条件的谱系进行 pseudotime 推断。
- `PseudotimeProjectionPlot()` 和 `DynamicHeatmap()`：展示 pseudotime 方向和谱系动态表达。

输出目录：

```text
results/GSE234129/figures/scop_extended/
results/GSE234129/tables/scop_extended/
results/GSE234129/reports/GSE234129_scop_extended_analysis_report.md
```

注意：

- 脚本固定 `OMP_NUM_THREADS`、`OPENBLAS_NUM_THREADS`、`MKL_NUM_THREADS`、`VECLIB_MAXIMUM_THREADS`、`BLIS_NUM_THREADS` 为 1，避免并行计算占满线程。
- `monocle3`、`tradeSeq`、`palantir` 当前未安装，正式轨迹分析使用已安装的 `slingshot`。
- `GroupHeatmap()` 在使用 `grouping.var` 时会按 `scop` 规则将表达比较限制为 `log2fc`。
- B/Plasma 的 Slingshot 轨迹可正常输出；该谱系动态热图可能因个别基因拟合出现 infinite/missing values 而被记录为局部跳过。
