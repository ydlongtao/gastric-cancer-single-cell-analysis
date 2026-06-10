# GSE234129 scDesign3 Simulation Report

Generated: 2026-06-10 09:51:05 EDT

## Inputs

- Annotated h5ad: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/objects/GSE234129_annotated.h5ad`
- Seurat counts: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/Datasets/GSE234129/processed/GSE234129_seurat.rds`
- QC baseline cells: 19144
- Feature strategy: HVG + classic markers + Leiden 0.5 top markers
- Feature count: 200
- max_features: 200

## scDesign3 Parameters

- seed: 20260609
- ncell: 76576
- celltype: `global_annotation`
- other_covariates: `sample`
- mu_formula: `global_annotation`
- sigma_formula: `1`
- family_use: `nb`
- usebam: `TRUE`
- copula: `gaussian`
- DT: `TRUE`

## Outputs

- `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_counts.rds`
- `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_sce.rds`
- `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_seurat.rds`
- `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/scdesign3/tables/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_metadata.tsv`
- `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/scdesign3/tables/GSE234129_scdesign3_original_vs_simulated_qc_summary.tsv`
- `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/scdesign3/tables/GSE234129_scdesign3_original_vs_simulated_composition.tsv`

## Validation

- Simulated cells: 76576
- Simulated genes: 200
- Metadata rows: 76576
- Status: passed

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
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ggplot2_4.0.3               scDesign3_1.10.0           
 [3] Seurat_5.5.0                SeuratObject_5.4.0         
 [5] sp_2.2-1                    SingleCellExperiment_1.34.0
 [7] SummarizedExperiment_1.42.0 Biobase_2.72.0             
 [9] GenomicRanges_1.64.0        Seqinfo_1.2.0              
[11] IRanges_2.46.0              S4Vectors_0.50.1           
[13] BiocGenerics_0.58.1         generics_0.1.4             
[15] MatrixGenerics_1.24.0       matrixStats_1.5.0          
[17] Matrix_1.7-5               

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3     jsonlite_2.0.0         magrittr_2.0.5        
  [4] spatstat.utils_3.2-3   farver_2.1.2           ragg_1.5.2            
  [7] vctrs_0.7.3            ROCR_1.0-12            spatstat.explore_3.8-1
 [10] htmltools_0.5.9        S4Arrays_1.12.0        SparseArray_1.12.2    
 [13] sctransform_0.4.3      parallelly_1.47.0      KernSmooth_2.23-26    
 [16] htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
 [19] plotly_4.12.0          zoo_1.8-15             igraph_2.3.2          
 [22] mime_0.13              lifecycle_1.0.5        pkgconfig_2.0.3       
 [25] R6_2.6.1               fastmap_1.2.0          fitdistrplus_1.2-6    
 [28] future_1.70.0          shiny_1.13.0           digest_0.6.39         
 [31] patchwork_1.3.2        tensor_1.5.1           RSpectra_0.16-2       
 [34] irlba_2.3.7            textshaping_1.0.5      labeling_0.4.3        
 [37] progressr_0.19.0       spatstat.sparse_3.2-0  httr_1.4.8            
 [40] polyclip_1.10-7        abind_1.4-8            mgcv_1.9-4            
 [43] compiler_4.6.0         withr_3.0.2            S7_0.2.2              
 [46] fastDummies_1.7.6      MASS_7.3-65            DelayedArray_0.38.2   
 [49] tools_4.6.0            lmtest_0.9-40          otel_0.2.0            
 [52] httpuv_1.6.17          future.apply_1.20.2    goftest_1.2-3         
 [55] glue_1.8.1             nlme_3.1-169           promises_1.5.0        
 [58] grid_4.6.0             Rtsne_0.17             cluster_2.1.8.2       
 [61] reshape2_1.4.5         gtable_0.3.6           spatstat.data_3.1-9   
 [64] tidyr_1.3.2            data.table_1.18.4      XVector_0.52.0        
 [67] spatstat.geom_3.8-1    RcppAnnoy_0.0.23       ggrepel_0.9.8         
 [70] sparseMVN_0.2.2        RANN_2.6.2             pillar_1.11.1         
 [73] stringr_1.6.0          spam_2.11-4            RcppHNSW_0.7.0        
 [76] later_1.4.8            splines_4.6.0          dplyr_1.2.1           
 [79] lattice_0.22-9         survival_3.8-6         deldir_2.0-4          
 [82] gamlss.data_6.0-7      tidyselect_1.2.1       miniUI_0.1.2          
 [85] pbapply_1.7-4          gridExtra_2.3          scattermore_1.2       
 [88] stringi_1.8.7          lazyeval_0.2.3         codetools_0.2-20      
 [91] tibble_3.3.1           cli_3.6.6              uwot_0.2.4            
 [94] systemfonts_1.3.2      xtable_1.8-8           reticulate_1.46.0     
 [97] Rcpp_1.1.1-1.1         globals_0.19.1         spatstat.random_3.5-0 
[100] png_0.1-9              spatstat.univar_3.2-0  parallel_4.6.0        
[103] dotCall64_1.2          mclust_6.1.2           listenv_0.10.1        
[106] gamlss.dist_6.1-1      viridisLite_0.4.3      mvtnorm_1.4-1         
[109] scales_1.4.0           ggridges_0.5.7         gamlss_5.5-0          
[112] purrr_1.2.2            rlang_1.2.0            cowplot_1.2.0         
```
