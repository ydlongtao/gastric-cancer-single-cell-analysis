from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


PROJECT_DIR = Path(__file__).resolve().parents[2]
INPUT_H5AD = PROJECT_DIR / "Datasets/GSE234129/processed/GSE234129_scanpy.h5ad"
OUT_BASE = PROJECT_DIR / "results/GSE234129"
QC_DIR = OUT_BASE / "qc"
FIG_DIR = OUT_BASE / "figures"
TABLE_DIR = OUT_BASE / "tables"
OBJECT_DIR = OUT_BASE / "objects"
OUT_H5AD = OBJECT_DIR / "GSE234129_qc_filtered.h5ad"

MIN_GENES = 300
MAX_GENES = 6000
MIN_COUNTS = 500
MAX_PCT_MT = 15


def ensure_dirs() -> None:
    for path in [QC_DIR, FIG_DIR, TABLE_DIR, OBJECT_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def save_counts_table(adata: sc.AnnData, prefix: str) -> None:
    for column in ["patient", "sample", "celltype"]:
        counts = adata.obs[column].value_counts().rename_axis(column).reset_index(name="n_cells")
        counts.to_csv(QC_DIR / f"{prefix}_{column}_counts.tsv", sep="\t", index=False)


def plot_qc(adata: sc.AnnData, prefix: str) -> None:
    metrics = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    for ax, metric in zip(axes, metrics, strict=True):
        values = adata.obs[metric].to_numpy()
        ax.hist(values, bins=80, color="#3a7ca5", edgecolor="white", linewidth=0.2)
        ax.set_title(metric)
        ax.set_ylabel("Cells")
    axes[0].axvline(MIN_COUNTS, color="#b84a39", linestyle="--", linewidth=1)
    axes[1].axvline(MIN_GENES, color="#b84a39", linestyle="--", linewidth=1)
    axes[1].axvline(MAX_GENES, color="#b84a39", linestyle="--", linewidth=1)
    axes[2].axvline(MAX_PCT_MT, color="#b84a39", linestyle="--", linewidth=1)
    fig.tight_layout()
    fig.savefig(QC_DIR / f"{prefix}_qc_histograms.png", dpi=180)
    plt.close(fig)


def summarize_qc(adata: sc.AnnData, prefix: str) -> pd.DataFrame:
    rows = []
    for metric in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        values = adata.obs[metric].to_numpy()
        rows.append(
            {
                "stage": prefix,
                "metric": metric,
                "min": float(np.min(values)),
                "q01": float(np.quantile(values, 0.01)),
                "q05": float(np.quantile(values, 0.05)),
                "median": float(np.median(values)),
                "mean": float(np.mean(values)),
                "q95": float(np.quantile(values, 0.95)),
                "q99": float(np.quantile(values, 0.99)),
                "max": float(np.max(values)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    adata = sc.read_h5ad(INPUT_H5AD)
    required = {"patient", "sample", "celltype"}
    missing = required.difference(adata.obs.columns)
    if missing:
        raise ValueError(f"Input object is missing required obs columns: {sorted(missing)}")

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    before_n = adata.n_obs
    save_counts_table(adata, "before_qc")
    plot_qc(adata, "before_qc")
    before_summary = summarize_qc(adata, "before_qc")

    keep = (
        (adata.obs["n_genes_by_counts"] >= MIN_GENES)
        & (adata.obs["n_genes_by_counts"] <= MAX_GENES)
        & (adata.obs["total_counts"] >= MIN_COUNTS)
        & (adata.obs["pct_counts_mt"] <= MAX_PCT_MT)
    )
    adata_qc = adata[keep].copy()

    after_n = adata_qc.n_obs
    save_counts_table(adata_qc, "after_qc")
    plot_qc(adata_qc, "after_qc")
    after_summary = summarize_qc(adata_qc, "after_qc")

    summary = pd.concat([before_summary, after_summary], ignore_index=True)
    summary.to_csv(QC_DIR / "GSE234129_qc_metric_summary.tsv", sep="\t", index=False)

    retention = []
    for column in ["patient", "sample", "celltype"]:
        before = adata.obs[column].value_counts()
        after = adata_qc.obs[column].value_counts()
        table = pd.DataFrame({"before_qc": before, "after_qc": after}).fillna(0).astype(int)
        table["retention_rate"] = table["after_qc"] / table["before_qc"]
        table.index.name = column
        table.reset_index().to_csv(QC_DIR / f"GSE234129_{column}_retention.tsv", sep="\t", index=False)
        retention.append((column, int((table["after_qc"] > 0).sum()), int(table.shape[0])))

    report_lines = [
        "# GSE234129 QC Report",
        "",
        f"Input: `{INPUT_H5AD}`",
        f"Output: `{OUT_H5AD}`",
        "",
        "## Thresholds",
        "",
        f"- n_genes_by_counts >= {MIN_GENES}",
        f"- n_genes_by_counts <= {MAX_GENES}",
        f"- total_counts >= {MIN_COUNTS}",
        f"- pct_counts_mt <= {MAX_PCT_MT}",
        "",
        "## Cell Counts",
        "",
        f"- Before QC: {before_n:,}",
        f"- After QC: {after_n:,}",
        f"- Removed: {before_n - after_n:,}",
        f"- Retention: {after_n / before_n:.2%}",
        "",
        "## Category Retention",
        "",
    ]
    report_lines.extend([f"- {name}: {kept}/{total} categories retained" for name, kept, total in retention])
    report_lines.append("")
    report_lines.append(f"Mitochondrial genes detected: {int(adata.var['mt'].sum())}")
    (QC_DIR / "GSE234129_qc_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    adata_qc.uns["qc_thresholds"] = {
        "min_genes": MIN_GENES,
        "max_genes": MAX_GENES,
        "min_counts": MIN_COUNTS,
        "max_pct_mt": MAX_PCT_MT,
    }
    adata_qc.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"Wrote {OUT_H5AD}")
    print(f"Cells before QC: {before_n}")
    print(f"Cells after QC: {after_n}")


if __name__ == "__main__":
    main()
