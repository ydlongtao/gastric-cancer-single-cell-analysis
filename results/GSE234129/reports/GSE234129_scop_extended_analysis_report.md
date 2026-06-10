# GSE234129 scop Extended Analysis Report

Generated: 2026-06-10 11:29:55 EDT

## Inputs

- QC scop object: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/objects/GSE234129_scop_plotting_seurat.rds`
- Raw scop object: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/objects/GSE234129_raw_scop_plotting_seurat.rds`
- Marker table: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/tables/GSE234129_leiden05_top20_markers.tsv`
- Cluster annotation table: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/tables/GSE234129_cluster_annotation.tsv`

## Object Summary

- QC dimensions: 27176 genes x 19144 cells
- Raw dimensions: 27176 genes x 19488 cells
- QC global annotations: 19
- QC Leiden 0.5 clusters: 20
- Raw Leiden 0.5 clusters: 19

## Completed Modules

- GroupHeatmap for classic markers, global annotation top markers, and Leiden top markers.
- FeatureHeatmap for sampled classic marker expression.
- CellCorHeatmap within the QC object and between QC/raw references.
- Slingshot trajectory inference for eligible lineage subsets.
- DynamicHeatmap attempted for each successful lineage with available marker/HVF features.

## Trajectory Dependency Status

- slingshot: `2.20.0`
- monocle3: `NOT_INSTALLED`
- tradeSeq: `NOT_INSTALLED`
- palantir: `NOT_INSTALLED`

## Outputs

- Figure directory: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/figures/scop_extended`
- Table directory: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/tables/scop_extended`
- Local trajectory object directory: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/objects/scop_extended`
- Step summary: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/tables/scop_extended/GSE234129_scop_extended_step_summary.tsv`
- Figure files written (PDF/TIF): 50
- Table files written: 10
- Local trajectory RDS objects written: 3

## Step Status

                                         step            status
                group_heatmap_classic_markers         completed
             group_heatmap_global_top_markers         completed
             group_heatmap_leiden_top_markers         completed
      feature_heatmap_classic_markers_sampled         completed
          cellcor_leiden_vs_global_annotation         completed
        cellcor_celltype_vs_global_annotation         completed
                    cellcor_sample_similarity         completed
 cellcor_raw_celltype_vs_qc_global_annotation         completed
              cellcor_raw_leiden_vs_qc_leiden         completed
              trajectory_t_nk_dynamic_heatmap         completed
                              trajectory_t_nk         completed
                      trajectory_wrapper_t_nk         completed
       trajectory_myeloid_tam_dynamic_heatmap         completed
                       trajectory_myeloid_tam         completed
               trajectory_wrapper_myeloid_tam         completed
          trajectory_b_plasma_dynamic_heatmap skipped_or_failed
                          trajectory_b_plasma         completed
                  trajectory_wrapper_b_plasma         completed
                  trajectory_epithelial_tumor           skipped
          trajectory_wrapper_epithelial_tumor         completed
            trajectory_endothelial_fibroblast           skipped
    trajectory_wrapper_endothelial_fibroblast         completed
                                                                                                                                detail
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
                                                                                                                                      
                                                                                           t_nk_Lineage1, t_nk_Lineage2, t_nk_Lineage3
                                                                                                                                      
                                                                                                                                      
                                                                      myeloid_tam_Lineage1, myeloid_tam_Lineage2, myeloid_tam_Lineage3
                                                                                                                                      
                                                                                    length of 'dimnames' [2] not equal to array extent
                                                                                                  b_plasma_Lineage1, b_plasma_Lineage2
                                                                                                                                      
                                                                                                  No cells matched the lineage pattern
                                                                                                                                      
 Fewer than 3 groups with at least 30 cells: Endothelial/Fibroblast | ref:Endothelial_C0=174; Endothelial/Fibroblast | ref:VSMC_C3=248
                                                                                                                                      

## Session Info

```text
R version 4.6.0 (2026-04-24)
Platform: aarch64-apple-darwin23
Running under: macOS Tahoe 26.5

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] future_1.70.0      patchwork_1.3.2    ggplot2_4.0.3      scop_0.8.9        
[5] Seurat_5.5.0       SeuratObject_5.4.0 sp_2.2-1          

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.23            TrajectoryUtils_1.20.0     
  [3] slingshot_2.20.0            splines_4.6.0              
  [5] later_1.4.8                 bitops_1.0-9               
  [7] tibble_3.3.1                polyclip_1.10-7            
  [9] fastDummies_1.7.6           lifecycle_1.0.5            
 [11] thisplot_0.4.0              doParallel_1.0.17          
 [13] globals_0.19.1              lattice_0.22-9             
 [15] MASS_7.3-65                 backports_1.5.1            
 [17] magrittr_2.0.5              plotly_4.12.0              
 [19] thisutils_0.4.7             httpuv_1.6.17              
 [21] otel_0.2.0                  sctransform_0.4.3          
 [23] spam_2.11-4                 spatstat.sparse_3.2-0      
 [25] reticulate_1.46.0           cowplot_1.2.0              
 [27] pbapply_1.7-4               RColorBrewer_1.1-3         
 [29] abind_1.4-8                 Rtsne_0.17                 
 [31] GenomicRanges_1.64.0        purrr_1.2.2                
 [33] BiocGenerics_0.58.1         tweenr_2.0.3               
 [35] circlize_0.4.18             IRanges_2.46.0             
 [37] S4Vectors_0.50.1            ggrepel_0.9.8              
 [39] irlba_2.3.7                 listenv_0.10.1             
 [41] spatstat.utils_3.2-3        goftest_1.2-3              
 [43] RSpectra_0.16-2             spatstat.random_3.5-0      
 [45] fitdistrplus_1.2-6          parallelly_1.47.0          
 [47] DelayedMatrixStats_1.34.0   codetools_0.2-20           
 [49] DelayedArray_0.38.2         RcppRoll_0.3.2             
 [51] ggforce_0.5.0               tidyselect_1.2.1           
 [53] shape_1.4.6.1               UCSC.utils_1.8.0           
 [55] farver_2.1.2                matrixStats_1.5.0          
 [57] stats4_4.6.0                spatstat.explore_3.8-1     
 [59] Seqinfo_1.2.0               jsonlite_2.0.0             
 [61] GetoptLong_1.1.1            progressr_0.19.0           
 [63] ggridges_0.5.7              survival_3.8-6             
 [65] iterators_1.0.14            systemfonts_1.3.2          
 [67] foreach_1.5.2               ggnewscale_0.5.2           
 [69] tools_4.6.0                 ragg_1.5.2                 
 [71] ica_1.0-3                   Rcpp_1.1.1-1.1             
 [73] glue_1.8.1                  gridExtra_2.3              
 [75] SparseArray_1.12.2          mgcv_1.9-4                 
 [77] MatrixGenerics_1.24.0       GenomeInfoDb_1.48.0        
 [79] dplyr_1.2.1                 withr_3.0.2                
 [81] fastmap_1.2.0               digest_0.6.39              
 [83] R6_2.6.1                    mime_0.13                  
 [85] textshaping_1.0.5           colorspace_2.1-2           
 [87] scattermore_1.2             Cairo_1.7-0                
 [89] tensor_1.5.1                spatstat.data_3.1-9        
 [91] tidyr_1.3.2                 generics_0.1.4             
 [93] data.table_1.18.4           httr_1.4.8                 
 [95] htmlwidgets_1.6.4           S4Arrays_1.12.0            
 [97] uwot_0.2.4                  pkgconfig_2.0.3            
 [99] gtable_0.3.6                ComplexHeatmap_2.28.0      
[101] lmtest_0.9-40               S7_0.2.2                   
[103] SingleCellExperiment_1.34.0 XVector_0.52.0             
[105] htmltools_0.5.9             dotCall64_1.2              
[107] clue_0.3-68                 scales_1.4.0               
[109] Biobase_2.72.0              png_0.1-9                  
[111] spatstat.univar_3.2-0       Signac_1.17.1              
[113] reshape2_1.4.5              rjson_0.2.23               
[115] checkmate_2.3.4             nlme_3.1-169               
[117] metR_0.18.3                 cachem_1.1.0               
[119] zoo_1.8-15                  GlobalOptions_0.1.4        
[121] stringr_1.6.0               KernSmooth_2.23-26         
[123] parallel_4.6.0              miniUI_0.1.2               
[125] pillar_1.11.1               grid_4.6.0                 
[127] proxyC_0.5.2                vctrs_0.7.3                
[129] RANN_2.6.2                  promises_1.5.0             
[131] xtable_1.8-8                cluster_2.1.8.2            
[133] princurve_2.1.6             magick_2.9.1               
[135] cli_3.6.6                   compiler_4.6.0             
[137] Rsamtools_2.28.0            rlang_1.2.0                
[139] crayon_1.5.3                future.apply_1.20.2        
[141] labeling_0.4.3              plyr_1.8.9                 
[143] stringi_1.8.7               viridisLite_0.4.3          
[145] deldir_2.0-4                BiocParallel_1.46.0        
[147] Biostrings_2.80.1           lazyeval_0.2.3             
[149] spatstat.geom_3.8-1         Matrix_1.7-5               
[151] RcppHNSW_0.7.0              sparseMatrixStats_1.24.0   
[153] shiny_1.13.0                SummarizedExperiment_1.42.0
[155] ROCR_1.0-12                 memoise_2.0.1              
[157] igraph_2.3.2                fastmatch_1.1-8            
```
