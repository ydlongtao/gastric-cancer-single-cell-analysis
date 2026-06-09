# GSE234129 Global Subcluster Annotation Report

Input: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/objects/GSE234129_clustered.h5ad`
Output: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/objects/GSE234129_annotated.h5ad`

Cluster key: `leiden_0.5`
Clusters annotated: 20

## Output Tables

- `GSE234129_leiden05_markers.tsv`
- `GSE234129_leiden05_celltype_crosstab.tsv`
- `GSE234129_classic_marker_scores.tsv`
- `GSE234129_cluster_annotation.tsv`

## Annotation Preview

|   cluster |   n_cells | dominant_existing_celltype   |   dominant_existing_celltype_fraction | classic_marker_based_label   | proposed_global_annotation                  | top10_marker_genes                                                         |
|----------:|----------:|:-----------------------------|--------------------------------------:|:-----------------------------|:--------------------------------------------|:---------------------------------------------------------------------------|
|         0 |      2995 | CD4_C1                       |                              0.370952 | T/NK                         | T/NK | ref:CD4_C1                           | LTB,EEF1A1,RPL32,RPL11,LDHB,RPL10,RPL13,RPL34,RPS12,RPS8                   |
|         1 |      2779 | CD8_C0                       |                              0.312343 | T/NK                         | T/NK | ref:CD8_C0                           | CCL5,CD8A,IL32,NKG7,GZMA,GZMK,CD8B,CD3E,PTPRCAP,CST7                       |
|        10 |       730 | TAM_C3                       |                              0.508219 | Myeloid/TAM                  | Myeloid/TAM | ref:TAM_C3                    | SAT1,PLAUR,FTH1,NAMPT,LYZ,BCL2A1,FCER1G,AIF1,SOD2,CST3                     |
|        11 |       544 | NK_C3                        |                              0.564338 | T/NK                         | T/NK | ref:NK_C3                            | CTSW,CD7,NKG7,XCL1,TRDC,IL2RB,KLRB1,XCL2,IFITM2,KLRD1                      |
|        12 |       383 | TAM_C4                       |                              0.584856 | Myeloid/TAM                  | Myeloid/TAM | ref:TAM_C4                    | HLA-DRA,CD74,HLA-DPB1,FTL,TMSB4X,HLA-DPA1,HLA-DQA1,C1QB,C1QC,C1QA          |
|        13 |       341 | Classical Mono_C1            |                              0.709677 | Myeloid/TAM                  | Myeloid/TAM | ref:Classical Mono_C1         | LYZ,S100A8,IER2,CTSS,S100A9,FCN1,S100A6,CCL3L1,CSTA,NFKBIA                 |
|        14 |       322 | CD8_C3                       |                              0.962733 | T/NK                         | T/NK | ref:CD8_C3                           | RPL32,RPS12,RPS21,RPL34,CD8B,RPL37,EEF1A1,RPS5,RPL19,RPS23                 |
|        15 |       266 | DNT                          |                              0.413534 | T/NK                         | T/NK | ref:DNT                              | ACTG1,GAPDH,ACTB,PPIA,PTMA,TUBB,CFL1,HMGB1,PFN1,STMN1                      |
|        16 |       248 | VSMC_C3                      |                              0.278226 | Endothelial/Fibroblast       | Endothelial/Fibroblast | ref:VSMC_C3        | CALD1,COL1A2,SPARC,COL3A1,IGFBP7,IFITM3,TPM2,COL6A2,COL1A1,TIMP1           |
|        17 |       174 | Endothelial_C0               |                              0.367816 | Endothelial/Fibroblast       | Endothelial/Fibroblast | ref:Endothelial_C0 | GNG11,IGFBP7,SPARC,IFI27,EGFL7,SPARCL1,IFITM3,IGFBP4,CRIP2,HSPG2           |
|        18 |       125 | Mast                         |                              0.984    | Myeloid/TAM                  | Myeloid/TAM | ref:Mast                      | TPSAB1,CPA3,TPSB2,MS4A2,CLU,KIT,VWA5A,GATA2,HPGDS,BACE2                    |
|        19 |        60 | pDCs                         |                              0.966667 | B/Plasma                     | B/Plasma | ref:pDCs                         | IRF7,SERPINF1,GZMB,IRF8,PLD4,UGCG,TCF4,TSPAN13,HERPUD1,SPIB                |
|         2 |      1788 | CD4_C2                       |                              0.291946 | T/NK                         | T/NK | ref:CD4_C2                           | TNFAIP3,JUNB,SPOCK2,RGS1,BTG1,DNAJB1,CD3D,ICOS,CD69,CREM                   |
|         3 |      1746 | CD8_C1                       |                              0.377434 | T/NK                         | T/NK | ref:CD8_C1                           | DUSP2,CD8A,NR4A2,CD69,CCL5,TNFAIP3,DNAJB1,CXCR4,CD3D,CCL4                  |
|         4 |      1429 | TAM_C0                       |                              0.843947 | Myeloid/TAM                  | Myeloid/TAM | ref:TAM_C0                    | CTSB,FTL,PSAP,CD68,CTSL,GLUL,C1QC,C1QA,CTSD,FCGRT                          |
|         5 |      1311 | Classical Mono_C1            |                              0.334859 | Myeloid/TAM                  | Myeloid/TAM | ref:Classical Mono_C1         | CST3,AIF1,TYROBP,LYZ,FTL,LST1,CD74,SERPINA1,HLA-DRA,HLA-DRB1               |
|         6 |      1279 | NK_C2                        |                              0.400313 | T/NK                         | T/NK | ref:NK_C2                            | NKG7,CST7,CCL4,GZMB,GNLY,NR4A2,PRF1,CTSW,CD247,FGFBP2                      |
|         7 |      1007 | Plasma_C0                    |                              0.48858  | B/Plasma                     | B/Plasma | ref:Plasma_C0                    | JCHAIN,SSR4,MZB1,DERL3,IGHA1,CD79A,SEC11C,FKBP11,HSP90B1,PRDX4             |
|         8 |       861 | B cell_C3                    |                              0.497096 | B/Plasma                     | B/Plasma | ref:B cell_C3                    | CD79A,MS4A1,CD37,CD74,HLA-DRA,HLA-DQA1,HLA-DRB1,HLA-DPB1,HLA-DQB1,HLA-DPA1 |
|         9 |       756 | NK_C1                        |                              0.535714 | T/NK                         | T/NK | ref:NK_C1                            | NKG7,PRF1,FGFBP2,GZMB,CTSW,GNLY,KLRD1,CST7,CD247,FCGR3A                    |
