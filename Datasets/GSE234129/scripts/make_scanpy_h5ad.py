from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


DATASET_DIR = Path(__file__).resolve().parents[1]
RAW_DIR = DATASET_DIR / "raw"
OUT_DIR = DATASET_DIR / "processed"


def main() -> None:
    matrix_path = RAW_DIR / "GSE234129_count_matrix.mtx.gz"
    barcodes_path = RAW_DIR / "GSE234129_barcodes.tsv.gz"
    features_path = RAW_DIR / "GSE234129_features.tsv.gz"
    meta_path = RAW_DIR / "GSE234129_meta.tsv.gz"
    out_path = OUT_DIR / "GSE234129_scanpy.h5ad"

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    counts = sc.read_mtx(matrix_path).T
    barcodes = pd.read_csv(barcodes_path, header=None, sep="\t")[0].astype(str)
    features = pd.read_csv(features_path, header=None, sep="\t")[0].astype(str)
    meta = pd.read_csv(meta_path, sep="\t")

    if counts.n_obs != len(barcodes):
        raise ValueError(f"Barcode count mismatch: matrix={counts.n_obs}, barcodes={len(barcodes)}")
    if counts.n_vars != len(features):
        raise ValueError(f"Feature count mismatch: matrix={counts.n_vars}, features={len(features)}")
    if "cell_barcodes" not in meta.columns:
        raise ValueError("Metadata file does not contain a 'cell_barcodes' column")

    counts.obs_names = barcodes.to_list()
    counts.var_names = features.to_list()
    counts.var_names_make_unique()

    meta = meta.astype({"cell_barcodes": str}).set_index("cell_barcodes")
    missing = counts.obs_names.difference(meta.index)
    if len(missing) > 0:
        raise ValueError(f"Metadata missing {len(missing)} cells; first missing cell: {missing[0]}")

    counts.obs = meta.loc[counts.obs_names].copy()
    counts.var["gene_symbol"] = counts.var_names
    counts.uns["source"] = "GEO GSE234129"
    counts.uns["raw_files"] = {
        "matrix": matrix_path.name,
        "barcodes": barcodes_path.name,
        "features": features_path.name,
        "metadata": meta_path.name,
    }

    ad.AnnData.write_h5ad(counts, out_path, compression="gzip")
    print(out_path)
    print(counts)


if __name__ == "__main__":
    main()
