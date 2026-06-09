from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import scanpy as sc


PROJECT_DIR = Path(__file__).resolve().parents[2]
IN_H5AD = PROJECT_DIR / "results/GSE234129/objects/GSE234129_qc_filtered.h5ad"
OUT_BASE = PROJECT_DIR / "results/GSE234129"
FIG_DIR = OUT_BASE / "figures"
TABLE_DIR = OUT_BASE / "tables"
OBJECT_DIR = OUT_BASE / "objects"
OUT_H5AD = OBJECT_DIR / "GSE234129_clustered.h5ad"

RESOLUTIONS = [0.3, 0.5, 0.8, 1.0]
MAIN_CLUSTER_KEY = "leiden_0.5"


def ensure_dirs() -> None:
    for path in [FIG_DIR, TABLE_DIR, OBJECT_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def main() -> None:
    ensure_dirs()
    adata = sc.read_h5ad(IN_H5AD)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=3000,
        flavor="seurat",
        batch_key="sample",
    )
    hvgs = int(adata.var["highly_variable"].sum())
    if hvgs == 0:
        raise ValueError("No highly variable genes were selected")

    sc.tl.pca(adata, n_comps=50, svd_solver="arpack", mask_var="highly_variable")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata, random_state=0)

    cluster_summary = []
    for resolution in RESOLUTIONS:
        key = f"leiden_{resolution:.1f}"
        sc.tl.leiden(adata, resolution=resolution, key_added=key, random_state=0)
        counts = adata.obs[key].value_counts().sort_index()
        for cluster, n_cells in counts.items():
            cluster_summary.append({"resolution": resolution, "cluster": cluster, "n_cells": int(n_cells)})

    pd.DataFrame(cluster_summary).to_csv(
        TABLE_DIR / "GSE234129_leiden_cluster_counts.tsv",
        sep="\t",
        index=False,
    )

    fig = sc.pl.umap(
        adata,
        color=MAIN_CLUSTER_KEY,
        legend_loc="on data",
        frameon=False,
        show=False,
        return_fig=True,
    )
    fig.savefig(FIG_DIR / "GSE234129_umap_leiden_0.5.png", dpi=220, bbox_inches="tight")

    for color in ["patient", "sample", "celltype"]:
        fig = sc.pl.umap(adata, color=color, frameon=False, show=False, return_fig=True)
        fig.savefig(FIG_DIR / f"GSE234129_umap_{color}.png", dpi=220, bbox_inches="tight")

    adata.uns["clustering_parameters"] = {
        "normalization_target_sum": 10000.0,
        "hvg_n_top_genes": 3000,
        "hvg_flavor": "seurat",
        "hvg_batch_key": "sample",
        "pca_n_comps": 50,
        "neighbors_n_neighbors": 15,
        "neighbors_n_pcs": 30,
        "leiden_resolutions": RESOLUTIONS,
        "main_cluster_key": MAIN_CLUSTER_KEY,
    }
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"Wrote {OUT_H5AD}")
    print(f"Highly variable genes: {hvgs}")
    print(pd.DataFrame(cluster_summary).groupby("resolution")["cluster"].nunique().to_string())


if __name__ == "__main__":
    main()
