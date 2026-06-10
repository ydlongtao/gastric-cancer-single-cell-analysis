from __future__ import annotations

from pathlib import Path

import scanpy as sc


PROJECT_DIR = Path(__file__).resolve().parents[2]
IN_H5AD = PROJECT_DIR / "Datasets/GSE234129/processed/GSE234129_scanpy.h5ad"
OUT_DIR = PROJECT_DIR / "results/GSE234129/tables"
OUT_TSV = OUT_DIR / "GSE234129_raw_scop_plot_metadata.tsv"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(IN_H5AD)
    required_obs = ["patient", "sample", "celltype"]
    missing = [column for column in required_obs if column not in adata.obs.columns]
    if missing:
        raise ValueError(f"Missing columns in raw h5ad: {missing}")

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat", batch_key="sample")
    sc.tl.pca(adata, n_comps=50, use_highly_variable=True, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata, random_state=20260610)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0.5", random_state=20260610)

    meta = adata.obs[
        [
            "patient",
            "sample",
            "celltype",
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "leiden_0.5",
        ]
    ].copy()
    meta.insert(0, "cell_barcode", adata.obs_names)
    meta["UMAP_1"] = adata.obsm["X_umap"][:, 0]
    meta["UMAP_2"] = adata.obsm["X_umap"][:, 1]
    meta.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"Wrote {OUT_TSV}")
    print(meta.shape)


if __name__ == "__main__":
    main()
