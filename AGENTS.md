# AGENTS.md

本文件用于指导后续在本项目中工作的自动化助手或协作者。

## 项目定位

这是一个胃癌单细胞转录组和空间转录组分析项目。默认目标是产出可复现的分析脚本、图表、表格和报告，而不是只给出操作建议。

## 主要语言和环境

- R：用于 Seurat、scDesign3、空间转录组统计和出版级图件。
- Python：用于 Scanpy、AnnData、数据格式转换、批处理和补充分析。
- Rust：用于后续可能需要的高性能辅助命令行工具；当前已安装 stable 工具链。

当前已确认的可用环境：

- Python/Scanpy 环境：`/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python`
- 系统 R：`/usr/local/bin/R` 和 `/usr/local/bin/Rscript`
- R 包：`Seurat` 5.5.0、`SeuratObject` 5.4.0、`scDesign3` 1.10.0
- Rust：`rustc` 1.96.0、`cargo` 1.96.0、`rustup` 1.29.0

默认 `/opt/homebrew/bin/python3` 未安装 `scanpy`。运行 Python 单细胞分析时，优先使用 `biomni_e1` 环境。

## 工作原则

- 先检查目录结构和已有文件，再创建或修改文件。
- 不覆盖原始数据；`data/raw/` 中的文件视为只读。
- 对大文件、原始测序数据、空间图像和中间对象保持谨慎，不随意复制或移动。
- 分析脚本应尽量可重复运行，避免依赖交互式 notebook 中的隐藏状态。
- 输出文件名应包含分析对象、方法或日期，避免使用 `final`、`new`、`test` 这类不可追踪名称。
- 修改已有分析时，保留旧结果或记录变更原因，除非用户明确要求清理。

## 建议目录约定

- `data/raw/`：原始数据。
- `data/processed/`：质控后或格式转换后的中间数据。
- `scripts/R/`：R 脚本。
- `scripts/python/`：Python 脚本。
- `scripts/rust/`：Rust 工具源码或实验。
- `notebooks/`：探索性 notebook。
- `results/figures/`：图件。
- `results/tables/`：表格。
- `results/reports/`：报告或阶段性总结。
- `docs/`：说明文档、数据字典、样本信息。

## 环境验证命令

检查 Scanpy：

```bash
/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python - <<'PY'
import scanpy as sc
print(sc.__version__)
PY
```

检查 Seurat 和 scDesign3：

```bash
Rscript -e 'library(Seurat); library(scDesign3); sessionInfo()'
```

检查 Rust：

```bash
cargo --version
rustc --version
cargo test --manifest-path scripts/rust/rust_check/Cargo.toml
```

## 代码风格

- R：优先使用清晰的函数封装；长流程脚本按质控、整合、注释、统计、绘图分段。
- Python：优先使用 `pathlib`、明确的输入输出路径和 `argparse`。
- Rust：如果加入 Rust 工具，优先做成小型 CLI，并在 README 中写明输入输出格式。
- 图件脚本应固定随机种子，并保存生成图件所需的关键参数。

## 图件输出规则

- 分析结果图件优先输出为出版和汇报友好的格式：`PDF` 和 `TIF/TIFF`。
- 同一张图建议同时保存一份矢量版 `.pdf` 和一份高分辨率位图 `.tif`。
- TIF 图默认使用至少 300 dpi；用于投稿或精细排版时优先使用 600 dpi。
- PNG 仅作为快速预览或网页展示格式，不作为主要交付格式。
- 图件应保存到 `results/.../figures/` 下的分析子目录，并使用能说明数据集、方法和分组的文件名。

## 单细胞和空间分析注意事项

- 明确记录物种、参考基因组、样本来源、分组信息和批次信息。
- 质控阈值不要硬编码在多个脚本中；建议集中记录在配置文件或脚本开头。
- 细胞类型注释需要保存 marker 依据、参考数据集或人工校正说明。
- 空间转录组分析需要记录坐标系统、切片编号、spot/cell 过滤标准和组织区域定义。
- 跨样本整合前先保留每个样本独立 QC 报告。

## 交付标准

完成一次分析或修改时，尽量给出：

- 变更的脚本或文档路径。
- 使用的环境和关键包版本。
- 生成的结果文件路径。
- 成功运行的验证命令。
- 未完成或需要用户确认的事项。
