from __future__ import annotations

from pathlib import Path

import scanpy as sc


PROJECT_DIR = Path(__file__).resolve().parents[2]
INPUT_H5AD = PROJECT_DIR / "results/GSE234129/objects/GSE234129_annotated.h5ad"
OUT_DIR = PROJECT_DIR / "results/GSE234129/scdesign3/tables"
METADATA_OUT = OUT_DIR / "GSE234129_scdesign3_input_metadata.tsv"
FEATURES_OUT = OUT_DIR / "GSE234129_scdesign3_input_features.tsv"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT_H5AD)

    metadata = adata.obs.copy()
    metadata.insert(0, "cell_barcode", adata.obs_names)
    metadata.to_csv(METADATA_OUT, sep="\t", index=False)

    features = adata.var.copy()
    features.insert(0, "gene", adata.var_names)
    features.to_csv(FEATURES_OUT, sep="\t", index=False)

    print(f"Wrote {METADATA_OUT}")
    print(f"Wrote {FEATURES_OUT}")
    print(f"Cells: {adata.n_obs}")
    print(f"Genes: {adata.n_vars}")


if __name__ == "__main__":
    main()
