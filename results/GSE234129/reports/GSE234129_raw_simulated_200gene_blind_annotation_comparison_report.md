# GSE234129 Raw + Simulated 200-Gene Blind Annotation Comparison

## Inputs

- Original Seurat object: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/Datasets/GSE234129/processed/GSE234129_seurat.rds`
- QC metadata: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/scdesign3/tables/GSE234129_scdesign3_input_metadata.tsv`
- Simulated Seurat object: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_seurat.rds`

## Method

- Merged 19144 QC-matched original cells and 76576 scDesign3 simulated cells on 200 shared genes.
- The clustering and blind annotation metadata contained only source/sample fields; original `celltype`, previous `global_annotation`, and simulated reference labels were reattached only after `blind_annotation` was assigned.
- SingleR/celldex HPCA completed at cluster level using genes shared with the 200-gene merged object.
- Marker panels were curated from the searched reference resources: CellMarker 2.0, PanglaoDB, TISCH2, and CellTypist documentation, with gastric cancer/TME-relevant broad lineages.

## Searched Annotation Resources

- SingleR/celldex HPCA: https://bioconductor.org/books/release/SingleRBook/introduction.html
- CellMarker 2.0: https://academic.oup.com/nar/article/51/D1/D870/6775381
- PanglaoDB: https://panglaodb.se/
- TISCH2: https://tisch.compbio.cn/documentation/
- CellTypist: https://celltypist.readthedocs.io/

## Cluster Annotation Summary

```text
 cluster n_cells    blind_annotation blind_annotation_confidence
       0   16783              T cell                        high
       1   14209              T cell                        high
      10    3241         Myeloid/TAM                        high
      11    2872             NK cell                        high
      12    1937         Myeloid/TAM                        high
      13    1223 Fibroblast/Pericyte                        high
      14     842              T cell                        high
      15     826         Myeloid/TAM                      medium
      16     797         Endothelial                        high
      17     350         Myeloid/TAM                        high
      18     321         Plasma cell                      medium
       2    9077              T cell                        high
       3    8809              T cell                        high
       4    7762         Myeloid/TAM                        high
       5    7600         Myeloid/TAM                        high
       6    5882             NK cell                        high
       7    5062         Plasma cell                      medium
       8    4292              B cell                        high
       9    3835             NK cell                        high
                         singler_hpca_label source_original_fraction
                 T_cell:CD4+_central_memory                0.1989513
                         T_cell:gamma-delta                0.1942431
                             Monocyte:CD16-                0.1940759
                                NK_cell:IL2                0.1960306
 Macrophage:monocyte-derived:IL-4/Dex/cntrl                0.2049561
                      iPS_cells:fibroblasts                0.2076860
                T_cell:CD4+_effector_memory                0.1199525
       Endothelial_cells:HUVEC:PR8-infected                0.3644068
       Endothelial_cells:lymphatic:TNFa_48h                0.1932246
                             Monocyte:CD14+                0.9628571
                                        GMP                0.1993769
                 T_cell:CD8+_Central_memory                0.1980831
                 T_cell:CD4+_central_memory                0.2028607
                             Monocyte:CD14+                0.1477712
 Macrophage:monocyte-derived:IL-4/Dex/cntrl                0.2130263
                                    NK_cell                0.2084325
                         B_cell:Plasma_cell                0.1981430
                            B_cell:immature                0.2029357
                                NK_cell:IL2                0.2070404
 source_simulated_fraction
                0.80104868
                0.80575691
                0.80592410
                0.80396936
                0.79504388
                0.79231398
                0.88004751
                0.63559322
                0.80677541
                0.03714286
                0.80062305
                0.80191693
                0.79713929
                0.85222881
                0.78697368
                0.79156749
                0.80185697
                0.79706431
                0.79295958
```

## Quantitative Comparison

```text
                                       metric     value
               original_blind_vs_celltype_ari 0.1195251
               original_blind_vs_celltype_nmi 0.5042635
   original_broad_blind_vs_celltype_broad_ari 0.8984599
   original_broad_blind_vs_celltype_broad_nmi 0.8280015
             simulated_blind_vs_reference_ari 0.3102584
             simulated_blind_vs_reference_nmi 0.6610355
 simulated_broad_blind_vs_reference_broad_ari 0.9989067
 simulated_broad_blind_vs_reference_broad_nmi 0.9961418
            mean_cluster_source_balance_score 0.1941231
          median_cluster_source_balance_score 0.1989513
```

## Interpretation

- Agreement should be interpreted primarily at broad lineage level because the selected scDesign3 simulation contains only 200 genes.
- Blind labels that collapse several held-out subtypes are expected when the old labels distinguish fine CD4/CD8/TAM/plasma states unsupported by the restricted gene set.
- Higher-confidence blind labels require both a clear marker-panel margin and agreement with SingleR HPCA where available.
- Low-confidence clusters should not be used for fine biological claims without rerunning scDesign3 with a larger gene set.

## Outputs

- Object: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/objects/GSE234129_raw_simulated_200gene_blind_annotated_seurat.rds`
- Tables: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/tables/blind_annotation`
- Figures: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/figures/blind_annotation`
