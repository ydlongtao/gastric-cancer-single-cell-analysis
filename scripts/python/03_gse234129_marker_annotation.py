from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc


PROJECT_DIR = Path(__file__).resolve().parents[2]
IN_H5AD = PROJECT_DIR / "results/GSE234129/objects/GSE234129_clustered.h5ad"
OUT_BASE = PROJECT_DIR / "results/GSE234129"
FIG_DIR = OUT_BASE / "figures"
TABLE_DIR = OUT_BASE / "tables"
OBJECT_DIR = OUT_BASE / "objects"
OUT_H5AD = OBJECT_DIR / "GSE234129_annotated.h5ad"

CLUSTER_KEY = "leiden_0.5"
ANNOTATION_KEY = "global_annotation"

MARKER_SETS = {
    "T/NK": ["CD3D", "CD3E", "NKG7", "GNLY", "CD8A", "IL7R"],
    "B/Plasma": ["MS4A1", "CD79A", "MZB1", "JCHAIN"],
    "Myeloid/TAM": ["LYZ", "LST1", "CST3", "CD68", "MRC1"],
    "Epithelial/Tumor-related": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "Endothelial/Fibroblast": ["PECAM1", "VWF", "COL1A1", "DCN"],
}


def ensure_dirs() -> None:
    for path in [FIG_DIR, TABLE_DIR, OBJECT_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def rank_genes_to_table(adata: sc.AnnData) -> pd.DataFrame:
    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)
    rows = []
    for group in groups:
        names = result["names"][group]
        scores = result["scores"][group]
        pvals = result["pvals"][group]
        pvals_adj = result["pvals_adj"][group]
        logfoldchanges = result["logfoldchanges"][group]
        for rank, gene in enumerate(names, start=1):
            rows.append(
                {
                    "cluster": group,
                    "rank": rank,
                    "gene": gene,
                    "score": float(scores[rank - 1]),
                    "logfoldchange": float(logfoldchanges[rank - 1]),
                    "pval": float(pvals[rank - 1]),
                    "pval_adj": float(pvals_adj[rank - 1]),
                }
            )
    return pd.DataFrame(rows)


def expression_matrix(adata: sc.AnnData, genes: list[str]) -> pd.DataFrame:
    present = [gene for gene in genes if gene in adata.var_names]
    if not present:
        return pd.DataFrame(index=adata.obs_names)
    matrix = adata[:, present].X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return pd.DataFrame(matrix, index=adata.obs_names, columns=present)


def marker_scores_by_cluster(adata: sc.AnnData) -> pd.DataFrame:
    rows = []
    for label, genes in MARKER_SETS.items():
        expr = expression_matrix(adata, genes)
        if expr.empty:
            scores = pd.Series(0.0, index=adata.obs_names)
            present = []
        else:
            scores = expr.mean(axis=1)
            present = list(expr.columns)
        tmp = pd.DataFrame({CLUSTER_KEY: adata.obs[CLUSTER_KEY].astype(str).to_numpy(), "score": scores.to_numpy()})
        cluster_scores = tmp.groupby(CLUSTER_KEY)["score"].mean()
        for cluster, score in cluster_scores.items():
            rows.append(
                {
                    "cluster": cluster,
                    "marker_set": label,
                    "mean_marker_score": float(score),
                    "genes_present": ",".join(present),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    ensure_dirs()
    adata = sc.read_h5ad(IN_H5AD)
    if CLUSTER_KEY not in adata.obs:
        raise ValueError(f"Missing cluster key: {CLUSTER_KEY}")

    sc.tl.rank_genes_groups(adata, groupby=CLUSTER_KEY, method="wilcoxon", use_raw=True)
    marker_table = rank_genes_to_table(adata)
    marker_table.to_csv(TABLE_DIR / "GSE234129_leiden05_markers.tsv", sep="\t", index=False)
    marker_table.query("rank <= 20").to_csv(
        TABLE_DIR / "GSE234129_leiden05_top20_markers.tsv",
        sep="\t",
        index=False,
    )

    cross = pd.crosstab(adata.obs[CLUSTER_KEY], adata.obs["celltype"])
    cross.to_csv(TABLE_DIR / "GSE234129_leiden05_celltype_crosstab.tsv", sep="\t")

    marker_scores = marker_scores_by_cluster(adata)
    marker_scores.to_csv(TABLE_DIR / "GSE234129_classic_marker_scores.tsv", sep="\t", index=False)

    annotations = []
    for cluster, row in cross.iterrows():
        cluster_total = int(row.sum())
        dominant_celltype = row.idxmax()
        dominant_fraction = float(row.max() / cluster_total)
        score_subset = marker_scores[marker_scores["cluster"] == str(cluster)].sort_values(
            "mean_marker_score", ascending=False
        )
        best_marker_set = score_subset.iloc[0]["marker_set"] if not score_subset.empty else "Unknown"
        top_genes = marker_table.loc[marker_table["cluster"] == str(cluster)].head(10)["gene"].tolist()
        annotations.append(
            {
                "cluster": cluster,
                "n_cells": cluster_total,
                "dominant_existing_celltype": dominant_celltype,
                "dominant_existing_celltype_fraction": dominant_fraction,
                "classic_marker_based_label": best_marker_set,
                "proposed_global_annotation": f"{best_marker_set} | ref:{dominant_celltype}",
                "top10_marker_genes": ",".join(top_genes),
            }
        )

    annotation_table = pd.DataFrame(annotations).sort_values("cluster")
    annotation_table.to_csv(TABLE_DIR / "GSE234129_cluster_annotation.tsv", sep="\t", index=False)

    mapping = dict(zip(annotation_table["cluster"].astype(str), annotation_table["proposed_global_annotation"], strict=True))
    adata.obs[ANNOTATION_KEY] = adata.obs[CLUSTER_KEY].astype(str).map(mapping).astype("category")

    marker_genes = sorted({gene for genes in MARKER_SETS.values() for gene in genes if gene in adata.var_names})
    if marker_genes:
        dotplot = sc.pl.dotplot(
            adata,
            marker_genes,
            groupby=CLUSTER_KEY,
            standard_scale="var",
            show=False,
            return_fig=True,
        )
        dotplot.savefig(FIG_DIR / "GSE234129_classic_marker_dotplot.png")

    fig = sc.pl.umap(adata, color=ANNOTATION_KEY, frameon=False, show=False, return_fig=True)
    fig.savefig(FIG_DIR / "GSE234129_umap_global_annotation.png", dpi=220, bbox_inches="tight")

    report_lines = [
        "# GSE234129 Global Subcluster Annotation Report",
        "",
        f"Input: `{IN_H5AD}`",
        f"Output: `{OUT_H5AD}`",
        "",
        f"Cluster key: `{CLUSTER_KEY}`",
        f"Clusters annotated: {annotation_table.shape[0]}",
        "",
        "## Output Tables",
        "",
        "- `GSE234129_leiden05_markers.tsv`",
        "- `GSE234129_leiden05_celltype_crosstab.tsv`",
        "- `GSE234129_classic_marker_scores.tsv`",
        "- `GSE234129_cluster_annotation.tsv`",
        "",
        "## Annotation Preview",
        "",
        annotation_table.to_markdown(index=False),
        "",
    ]
    (OUT_BASE / "GSE234129_analysis_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    adata.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"Wrote {OUT_H5AD}")
    print(f"Wrote {TABLE_DIR / 'GSE234129_cluster_annotation.tsv'}")
    print(f"Clusters annotated: {annotation_table.shape[0]}")


if __name__ == "__main__":
    main()
