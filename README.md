# Gastric Cancer Single-Cell Analysis

中文说明：胃癌单细胞项目

本项目用于整合分析胃癌相关的单细胞转录组和空间转录组数据。主要语言和工具链包括 R、Python、Rust 及其生态程序包。

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
