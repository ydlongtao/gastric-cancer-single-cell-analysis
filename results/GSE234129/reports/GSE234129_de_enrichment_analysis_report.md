# GSE234129 Differential Expression and Enrichment Report

Generated: 2026-06-10 11:58:20 EDT

## Inputs

- Marker table: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/tables/GSE234129_leiden05_markers.tsv`
- Cluster annotation: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/tables/GSE234129_cluster_annotation.tsv`
- Seurat plotting object: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/objects/GSE234129_scop_plotting_seurat.rds`

## Differential Expression

- Source DE statistics are the existing Scanpy Wilcoxon Leiden 0.5 marker results.
- Significant markers are defined as `pval_adj < 0.05` and `abs(logfoldchange) >= 0.25`.
- Technical genes matching mitochondrial, ribosomal, or hemoglobin prefixes are excluded from enrichment gene sets.
- Tested marker rows: 543520
- Unique tested genes: 27176
- Mapped Entrez background genes: 19364
- Significant up marker rows used for ORA: 26829

## Enrichment

- GO Biological Process enrichment was run per `global_annotation` and per broad annotation.
- KEGG, Reactome, and MSigDB Hallmark ORA were run per broad annotation.
- ORA background is all tested genes from the marker table with Entrez mappings.

## Outputs

- Figure directory: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/figures/de_enrichment`
- Table directory: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/tables/de_enrichment`
- Step summary: `/Users/huangfulongtao/Desktop/胃癌单细胞项目/results/GSE234129/tables/de_enrichment/GSE234129_de_enrichment_step_summary.tsv`
- Figure files written (PDF/TIF): 16
- Table files written: 12

## Step Status

                             step    status detail
         volcano_faceted_clusters completed       
               de_summary_barplot completed       
                de_marker_heatmap completed       
   enrich_go_bp_global_annotation completed       
     enrich_kegg_broad_annotation completed       
 enrich_reactome_broad_annotation completed       
 enrich_hallmark_broad_annotation completed       
   enrichment_network_go_bp_broad completed       

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
 [1] scop_0.8.9             msigdbr_26.1.0         ReactomePA_1.56.0     
 [4] enrichplot_1.32.0      org.Hs.eg.db_3.23.1    AnnotationDbi_1.74.0  
 [7] IRanges_2.46.0         S4Vectors_0.50.1       Biobase_2.72.0        
[10] BiocGenerics_0.58.1    generics_0.1.4         clusterProfiler_4.20.0
[13] ggrepel_0.9.8          ggplot2_4.0.3          Seurat_5.5.0          
[16] SeuratObject_5.4.0     sp_2.2-1              

loaded via a namespace (and not attached):
  [1] bitops_1.0-9             fs_2.1.0                 matrixStats_1.5.0       
  [4] spatstat.sparse_3.2-0    httr_1.4.8               RColorBrewer_1.1-3      
  [7] doParallel_1.0.17        tools_4.6.0              sctransform_0.4.3       
 [10] R6_2.6.1                 lazyeval_0.2.3           uwot_0.2.4              
 [13] GetoptLong_1.1.1         withr_3.0.2              graphite_1.58.0         
 [16] gridExtra_2.3            progressr_0.19.0         textshaping_1.0.5       
 [19] cli_3.6.6                Cairo_1.7-0              spatstat.explore_3.8-1  
 [22] fastDummies_1.7.6        scatterpie_0.2.6         labeling_0.4.3          
 [25] S7_0.2.2                 spatstat.data_3.1-9      ggridges_0.5.7          
 [28] pbapply_1.7-4            Rsamtools_2.28.0         systemfonts_1.3.2       
 [31] yulab.utils_0.2.4        gson_0.1.0               DOSE_4.6.0              
 [34] parallelly_1.47.0        RSQLite_3.53.1           gridGraphics_0.5-1      
 [37] shape_1.4.6.1            ica_1.0-3                spatstat.random_3.5-0   
 [40] dplyr_1.2.1              GO.db_3.23.1             Matrix_1.7-5            
 [43] abind_1.4-8              lifecycle_1.0.5          yaml_2.3.12             
 [46] qvalue_2.44.0            Rtsne_0.17               grid_4.6.0              
 [49] blob_1.3.0               thisplot_0.4.0           promises_1.5.0          
 [52] crayon_1.5.3             miniUI_0.1.2             ggtangle_0.1.2          
 [55] lattice_0.22-9           cowplot_1.2.0            KEGGREST_1.52.0         
 [58] magick_2.9.1             pillar_1.11.1            ComplexHeatmap_2.28.0   
 [61] GenomicRanges_1.64.0     rjson_0.2.23             future.apply_1.20.2     
 [64] codetools_0.2-20         fastmatch_1.1-8          glue_1.8.1              
 [67] ggiraph_0.9.6            ggfun_0.2.0              spatstat.univar_3.2-0   
 [70] fontLiberation_0.1.0     data.table_1.18.4        vctrs_0.7.3             
 [73] png_0.1-9                treeio_1.36.1            spam_2.11-4             
 [76] gtable_0.3.6             assertthat_0.2.1         cachem_1.1.0            
 [79] Signac_1.17.1            mime_0.13                tidygraph_1.3.1         
 [82] Seqinfo_1.2.0            survival_3.8-6           aisdk_1.4.12            
 [85] RcppRoll_0.3.2           iterators_1.0.14         thisutils_0.4.7         
 [88] fitdistrplus_1.2-6       ROCR_1.0-12              nlme_3.1-169            
 [91] ggtree_4.2.0             bit64_4.8.2              fontquiver_0.2.1        
 [94] RcppAnnoy_0.0.23         GenomeInfoDb_1.48.0      irlba_2.3.7             
 [97] KernSmooth_2.23-26       otel_0.2.0               colorspace_2.1-2        
[100] DBI_1.3.0                tidyselect_1.2.1         processx_3.9.0          
[103] bit_4.6.0                compiler_4.6.0           curl_7.1.0              
[106] httr2_1.2.2              graph_1.90.0             fontBitstreamVera_0.1.1 
[109] plotly_4.12.0            scales_1.4.0             lmtest_0.9-40           
[112] callr_3.8.0              rappdirs_0.3.4           stringr_1.6.0           
[115] digest_0.6.39            goftest_1.2-3            spatstat.utils_3.2-3    
[118] XVector_0.52.0           htmltools_0.5.9          pkgconfig_2.0.3         
[121] MatrixGenerics_1.24.0    sparseMatrixStats_1.24.0 fastmap_1.2.0           
[124] UCSC.utils_1.8.0         rlang_1.2.0              GlobalOptions_0.1.4     
[127] htmlwidgets_1.6.4        shiny_1.13.0             farver_2.1.2            
[130] zoo_1.8-15               jsonlite_2.0.0           BiocParallel_1.46.0     
[133] GOSemSim_2.38.0          magrittr_2.0.5           ggplotify_0.1.3         
[136] dotCall64_1.2            patchwork_1.3.2          Rcpp_1.1.1-1.1          
[139] ape_5.8-1                babelgene_22.9           ggnewscale_0.5.2        
[142] viridis_0.6.5            gdtools_0.5.1            reticulate_1.46.0       
[145] stringi_1.8.7            ggraph_2.2.2             MASS_7.3-65             
[148] plyr_1.8.9               parallel_4.6.0           listenv_0.10.1          
[151] deldir_2.0-4             Biostrings_2.80.1        graphlayouts_1.2.3      
[154] splines_4.6.0            tensor_1.5.1             circlize_0.4.18         
[157] igraph_2.3.2             spatstat.geom_3.8-1      enrichit_0.1.4          
[160] RcppHNSW_0.7.0           reshape2_1.4.5           foreach_1.5.2           
[163] tweenr_2.0.3             httpuv_1.6.17            RANN_2.6.2              
[166] tidyr_1.3.2              purrr_1.2.2              polyclip_1.10-7         
[169] future_1.70.0            clue_0.3-68              scattermore_1.2         
[172] ggforce_0.5.0            xtable_1.8-8             reactome.db_1.96.0      
[175] RSpectra_0.16-2          tidytree_0.4.7           tidydr_0.0.6            
[178] later_1.4.8              ragg_1.5.2               viridisLite_0.4.3       
[181] tibble_3.3.1             aplot_0.2.9              memoise_2.0.1           
[184] cluster_2.1.8.2          globals_0.19.1          
```
