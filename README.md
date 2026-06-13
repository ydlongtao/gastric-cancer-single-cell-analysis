# Gastric Cancer Single-Cell Analysis

中文说明：胃癌单细胞项目

本项目用于整合分析胃癌相关的单细胞转录组和空间转录组数据。主要语言和工具链包括 R、Python、Rust 及其生态程序包。

## 项目亮点

本仓库围绕胃癌单细胞数据集 `GSE234129` 建立了一套可复现的分析流程，覆盖原始数据整理、Scanpy/Seurat 质控聚类、scop 可视化、差异表达和通路富集、scDesign3 模拟，以及模拟数据与原始数据的系统比较。

当前重点结果是：基于 19,144 个 QC 后原始细胞，用 `scDesign3` 生成 76,576 个模拟细胞，并在 200 个共同基因上将原始数据与模拟数据合并，进行不参考原有注释的盲法细胞群重新注释。结果显示，200 基因合并对象能够较好恢复 broad lineage 层面的结构，但不适合直接作为 CD4/CD8/TAM/plasma 等精细亚群的最终判定。

关键比较指标：

| 比较对象 | ARI | NMI | 解读 |
| --- | ---: | ---: | --- |
| 原始数据：盲法注释 vs 原有精细 celltype | 0.120 | 0.504 | 精细亚群一致性有限，符合 200 基因输入的预期限制。 |
| 原始数据：盲法 broad lineage vs 原有 broad lineage | 0.898 | 0.828 | 大类谱系恢复较好。 |
| 模拟数据：盲法注释 vs 参考标签 | 0.310 | 0.661 | 精细标签层面有合并和塌缩。 |
| 模拟数据：盲法 broad lineage vs 参考 broad lineage | 0.999 | 0.996 | 模拟数据在大类谱系层面高度一致。 |

## 图件总览

### 1. 主流程细胞注释和整体结构

QC 后细胞在 UMAP 空间中形成 T/NK、B/Plasma、Myeloid/TAM、Stromal 等主要结构。

<p align="center">
  <img src="docs/readme_figures/gse234129_scop_umap_global_annotation.png" alt="GSE234129 global annotation UMAP" width="760">
</p>

### 2. scDesign3 原始数据 vs 模拟数据

scDesign3 四倍模拟用于考察当前分析流程在生成数据上的稳定性。UMAP 显示原始细胞和模拟细胞整体能够投射到相近表达空间；组成和 marker 表达图用于进一步检查模拟结果是否保留主要细胞群结构。

| 原始/模拟 UMAP | 细胞组成对比 |
| --- | --- |
| <img src="docs/readme_figures/gse234129_scdesign3_original_vs_simulated_umap.png" alt="Original vs simulated UMAP" width="420"> | <img src="docs/readme_figures/gse234129_scdesign3_original_vs_simulated_composition.png" alt="Original vs simulated composition" width="420"> |

| 原始/模拟 marker 表达 |
| --- |
| <img src="docs/readme_figures/gse234129_scdesign3_original_vs_simulated_marker_dotplot.png" alt="Original vs simulated marker dotplot" width="860"> |

### 3. 原始数据与模拟数据合并后的盲法重新注释

为避免直接继承原有注释，本项目将原始细胞与模拟细胞合并后重新聚类，并使用 marker 面板和 `SingleR`/`celldex` HPCA 进行 blind annotation。原有 `celltype` 和模拟参考标签仅在最后作为 holdout label 参与比较。

| 合并对象按来源着色 | 合并对象按盲法注释着色 |
| --- | --- |
| <img src="docs/readme_figures/gse234129_blind_umap_source.png" alt="Merged raw and simulated UMAP by source" width="420"> | <img src="docs/readme_figures/gse234129_blind_umap_annotation.png" alt="Merged raw and simulated UMAP by blind annotation" width="420"> |

| 原始细胞：盲法注释 vs 原有 celltype | 模拟细胞：盲法注释 vs 参考标签 |
| --- | --- |
| <img src="docs/readme_figures/gse234129_blind_original_heatmap.png" alt="Blind annotation versus original celltype" width="420"> | <img src="docs/readme_figures/gse234129_blind_simulated_heatmap.png" alt="Blind annotation versus simulated reference labels" width="420"> |

| 盲法注释 marker dotplot |
| --- |
| <img src="docs/readme_figures/gse234129_blind_marker_dotplot.png" alt="Blind annotation marker dotplot" width="860"> |

重点解释：

- `Cluster 15`：marker 面板偏向 `Myeloid/TAM`，但 `SingleR` 指向 endothelial，属于本次最需要谨慎解释的混合/低信息 cluster。
- `Cluster 17`：由 `S100A8`、`FCN1`、`S100A9`、`LYZ` 等支持，且原始细胞占比 96.3%，提示可能是真实数据中特异的炎症性单核/髓系状态。
- `Cluster 18`：`JCHAIN`、`MZB1`、`DERL3` 支持 plasma cell 倾向，但 `SingleR` 映射异常，不建议做精细亚型判断。

相关报告：

- [中文 HTML 解读报告](results/GSE234129/reports/GSE234129_raw_simulated_200gene_blind_annotation_chinese_report.html)
- [English HTML report](results/GSE234129/reports/GSE234129_raw_simulated_200gene_blind_annotation_english_report.html)
- [Markdown analysis report](results/GSE234129/reports/GSE234129_raw_simulated_200gene_blind_annotation_comparison_report.md)

### 4. 差异表达和通路富集

在 Leiden 0.5 cluster 和全局注释基础上整理差异表达 marker，并进行 GO、KEGG、Reactome 和 Hallmark ORA 富集。下图展示不同全局注释群体的 GO Biological Process 富集结果。

<p align="center">
  <img src="docs/readme_figures/gse234129_de_go_bp_dotplot_global_annotation.png" alt="GO BP enrichment dotplot by global annotation" width="760">
</p>

## 研究内容

- 单细胞转录组数据质控、标准化、降维、聚类和细胞类型注释。
- 空间转录组数据质控、空间域识别、组织结构关联和可视化。
- 单细胞与空间转录组联合分析，包括细胞类型映射、空间去卷积和肿瘤微环境解析。
- 使用模拟、重采样或生成模型评估分析流程稳定性，例如 `scDesign3`。
- 在需要性能优化或可复用命令行工具时，使用 Rust 编写辅助程序。

## 当前环境检查

检查日期：2026-06-09

### Python

默认 Python：

```bash
/opt/homebrew/bin/python3 --version
# Python 3.14.5
```

默认 Python 环境中暂未发现 `scanpy`、`anndata`、`numpy`、`pandas`。

已发现可用的 micromamba 环境：

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python
```

该环境中已部署：

- `scanpy` 1.11.5
- `anndata` 0.12.16
- `numpy` 2.1.0
- `pandas` 2.3.3

建议使用方式：

```bash
/Users/huangfulongtao/.local/bin/micromamba activate biomni_e1
python -c "import scanpy as sc; print(sc.__version__)"
```

如果当前 shell 尚未初始化 micromamba，也可以直接使用环境内解释器：

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/example_scanpy.py
```

### R

系统 R：

```bash
/usr/local/bin/R --version
# R version 4.6.0 (2026-04-24)
```

系统 R 中已部署：

- `Seurat` 5.5.0
- `SeuratObject` 5.4.0
- `scDesign3` 1.10.0

验证命令：

```bash
Rscript -e 'pkgs <- c("Seurat", "scDesign3"); for (p in pkgs) print(packageVersion(p))'
```

### Rust

已通过 Homebrew 安装 `rustup`，并安装 stable Rust 工具链。

已在 `~/.zshrc` 中加入：

```bash
export PATH="/opt/homebrew/opt/rustup/bin:$PATH"
```

当前版本：

- `rustup` 1.29.0
- `rustc` 1.96.0
- `cargo` 1.96.0
- `rustfmt` 1.9.0
- `clippy` 0.1.96

验证项目：

```bash
scripts/rust/rust_check/
```

验证命令：

```bash
cargo test --manifest-path scripts/rust/rust_check/Cargo.toml
cargo run --manifest-path scripts/rust/rust_check/Cargo.toml
```

## 建议目录结构

```text
.
├── README.md
├── AGENTS.md
├── data/
│   ├── raw/
│   ├── processed/
│   └── reference/
├── scripts/
│   ├── R/
│   ├── python/
│   └── rust/
├── notebooks/
├── results/
│   ├── figures/
│   ├── tables/
│   └── reports/
└── docs/
```

## 数据管理约定

- `Datasets/` 保存项目使用的公共数据集、单细胞/空间转录组对象和病历资料。
- `data/raw/` 保存原始数据，原则上只读，不直接覆盖。
- `data/processed/` 保存中间对象，例如 `.h5ad`、`.rds`、`.h5seurat`、`.loom`。
- `results/figures/` 保存图件，文件名建议包含分析步骤、数据集和日期。
- `results/tables/` 保存差异基因、marker、细胞比例、空间域统计等表格。
- 大体积数据不建议提交到版本库；必要时使用清单文件记录来源、下载链接、校验值和处理脚本。

当前已整理数据集：

- `Datasets/GSE234129/processed/GSE234129_scanpy.h5ad`
- `Datasets/GSE234129/processed/GSE234129_seurat.rds`

当前已完成 Scanpy 主流程分析：

- `results/GSE234129/objects/GSE234129_qc_filtered.h5ad`
- `results/GSE234129/objects/GSE234129_clustered.h5ad`
- `results/GSE234129/objects/GSE234129_annotated.h5ad`
- `results/GSE234129/GSE234129_analysis_report.md`

## 常用启动命令

Python/Scanpy：

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python -c "import scanpy as sc; print(sc.__version__)"
```

GSE234129 Scanpy 分析流程：

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/01_gse234129_qc.py
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/02_gse234129_clustering.py
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/03_gse234129_marker_annotation.py
```

GSE234129 scop 绘图流程：

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/export_gse234129_scop_inputs.py
Rscript scripts/R/04_gse234129_scop_plots.R
```

GSE234129 原始未过滤数据集 scop 绘图流程：

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python scripts/python/export_gse234129_raw_scop_inputs.py
Rscript scripts/R/06_gse234129_raw_scop_plots.R
```

原始数据集 scop 输出：

- 输入细胞数：19,488
- sample 数：17
- 原始 celltype 标签数：62
- Leiden 0.5 cluster 数：19
- 图件目录：`results/GSE234129/figures/scop_raw/`
- 正式图件格式：PDF 和 600 dpi TIF；PNG 仅作预览

GSE234129 scop 扩展分析流程：

```bash
Rscript scripts/R/08_gse234129_scop_extended_analysis.R
```

扩展分析内容：

- `GroupHeatmap`：经典 marker、`global_annotation` top markers、Leiden 0.5 top markers。
- `CellCorHeatmap`：QC 对象内部的 cluster/annotation/sample 相似性，以及 raw-vs-QC 对照。
- Slingshot trajectory inference：T/NK、Myeloid/TAM、B/Plasma 三条满足规模条件的谱系。
- `PseudotimeProjectionPlot` 和动态基因热图：用于展示 pseudotime 方向和谱系相关表达变化。

扩展分析输出：

- 图件目录：`results/GSE234129/figures/scop_extended/`
- 表格目录：`results/GSE234129/tables/scop_extended/`
- 报告：`results/GSE234129/reports/GSE234129_scop_extended_analysis_report.md`
- 正式图件格式：PDF 和 600 dpi TIF；PNG 仅作预览

GSE234129 差异表达和富集分析流程：

```bash
Rscript scripts/R/09_gse234129_de_enrichment_analysis.R
```

分析内容：

- 基于既有 Scanpy Wilcoxon Leiden 0.5 marker 结果整理差异表达表。
- 输出 cluster 分面火山图、显著 marker 数量图、top marker group heatmap。
- 对显著上调 marker 进行 GO Biological Process、KEGG、Reactome 和 MSigDB Hallmark ORA 富集。
- 使用所有已测试且可映射 Entrez ID 的 marker 基因作为 ORA background。

输出：

- 图件目录：`results/GSE234129/figures/de_enrichment/`
- 表格目录：`results/GSE234129/tables/de_enrichment/`
- 报告：`results/GSE234129/reports/GSE234129_de_enrichment_analysis_report.md`
- 完整注释 marker 表较大，保留在本地；Git 中同步摘要、top marker、富集表、报告和 PDF/TIF 图件。

GSE234129 scDesign3 四倍模拟流程：

```bash
Rscript scripts/R/05_gse234129_scdesign3_simulation.R \
  --refresh-inputs \
  --max-features=200 \
  --n-cores=7 \
  --mu-formula=global_annotation
```

当前成功模拟版本：

- 输入对象：`results/GSE234129/objects/GSE234129_annotated.h5ad`
- 基准细胞数：19,144
- 模拟细胞数：76,576
- 特征基因数：200
- 模型公式：`global_annotation`
- 运行耗时：25 分 48 秒（2026-06-10 09:25:17 到 09:51:05）
- 报告：`results/GSE234129/scdesign3/GSE234129_scdesign3_simulation_report.md`

R/Seurat：

```bash
Rscript -e 'library(Seurat); sessionInfo()'
```

R/scDesign3：

```bash
Rscript -e 'library(scDesign3); packageVersion("scDesign3")'
```

Rust/Cargo：

```bash
cargo --version
cargo run --manifest-path scripts/rust/rust_check/Cargo.toml
```

## 可复现性记录

每个分析脚本建议在输出目录中保存：

- 输入数据路径和数据版本。
- 关键参数。
- 使用的软件版本。
- 生成时间。
- 输出文件列表。
