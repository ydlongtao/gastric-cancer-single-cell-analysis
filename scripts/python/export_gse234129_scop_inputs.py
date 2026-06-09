from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc


PROJECT_DIR = Path(__file__).resolve().parents[2]
IN_H5AD = PROJECT_DIR / "results/GSE234129/objects/GSE234129_annotated.h5ad"
OUT_DIR = PROJECT_DIR / "results/GSE234129/tables"
OUT_TSV = OUT_DIR / "GSE234129_scop_plot_metadata.tsv"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(IN_H5AD)
    required_obs = [
        "patient",
        "sample",
        "celltype",
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "leiden_0.5",
        "global_annotation",
    ]
    missing = [column for column in required_obs if column not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing columns in annotated h5ad: {missing}")
    if "X_umap" not in adata.obsm:
        raise ValueError("Missing UMAP coordinates in adata.obsm['X_umap']")

    meta = adata.obs[required_obs].copy()
    meta.insert(0, "cell_barcode", adata.obs_names)
    meta["UMAP_1"] = adata.obsm["X_umap"][:, 0]
    meta["UMAP_2"] = adata.obsm["X_umap"][:, 1]
    meta.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Wrote {OUT_TSV}")
    print(meta.shape)


if __name__ == "__main__":
    main()
