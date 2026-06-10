suppressPackageStartupMessages({
  library(Seurat)
  library(scop)
  library(ggplot2)
  library(patchwork)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

script_arg <- commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))]
if (length(script_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
  project_dir <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
} else {
  project_dir <- getwd()
}

message_time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

qc_object_path <- file.path(project_dir, "results/GSE234129/objects/GSE234129_scop_plotting_seurat.rds")
raw_object_path <- file.path(project_dir, "results/GSE234129/objects/GSE234129_raw_scop_plotting_seurat.rds")
marker_path <- file.path(project_dir, "results/GSE234129/tables/GSE234129_leiden05_top20_markers.tsv")
cluster_annotation_path <- file.path(project_dir, "results/GSE234129/tables/GSE234129_cluster_annotation.tsv")

figure_dir <- file.path(project_dir, "results/GSE234129/figures/scop_extended")
table_dir <- file.path(project_dir, "results/GSE234129/tables/scop_extended")
report_dir <- file.path(project_dir, "results/GSE234129/reports")
object_dir <- file.path(project_dir, "results/GSE234129/objects/scop_extended")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(table_dir, "GSE234129_scop_extended_step_summary.tsv")
report_path <- file.path(report_dir, "GSE234129_scop_extended_analysis_report.md")

classic_markers <- c(
  "CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CCR7", "LTB",
  "CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "NKG7", "GNLY", "PRF1",
  "MS4A1", "CD79A", "CD79B", "CD74", "MZB1", "JCHAIN", "IGHG1",
  "LYZ", "LST1", "S100A8", "S100A9", "CST3", "AIF1", "FCGR3A",
  "CD68", "C1QA", "C1QB", "C1QC", "MRC1", "FCER1G",
  "EPCAM", "KRT8", "KRT18", "KRT19", "MUC1",
  "PECAM1", "VWF", "ENG", "KDR",
  "COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "RGS5",
  "MKI67", "TOP2A", "STMN1", "TPSAB1", "CPA3", "IRF7", "GZMB"
)

technical_gene <- function(x) {
  grepl("^MT-|^RPS|^RPL|^MTRNR|^HB[ABDEGQMZ]", x)
}

read_tsv <- function(path) {
  if (!file.exists(path)) {
    stop("Missing input file: ", path)
  }
  read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  message_time("Wrote ", path)
}

step_results <- data.frame(
  step = character(),
  status = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

record_status <- function(step, status, detail = "") {
  step_results <<- rbind(
    step_results,
    data.frame(step = step, status = status, detail = detail, stringsAsFactors = FALSE)
  )
}

run_step <- function(step, expr) {
  tryCatch(
    {
      force(expr)
      record_status(step, "completed", "")
    },
    error = function(e) {
      msg <- conditionMessage(e)
      record_status(step, "skipped_or_failed", msg)
      message_time("Skipped/failed ", step, ": ", msg)
    }
  )
}

save_plot <- function(plot, basename, width = 8, height = 6, dpi = 600) {
  pdf_path <- file.path(figure_dir, paste0(basename, ".pdf"))
  tif_path <- file.path(figure_dir, paste0(basename, ".tif"))
  png_path <- file.path(figure_dir, paste0(basename, ".png"))
  ggsave(pdf_path, plot = plot, width = width, height = height, device = cairo_pdf)
  ggsave(tif_path, plot = plot, width = width, height = height, dpi = dpi, device = "tiff", compression = "lzw")
  ggsave(png_path, plot = plot, width = width, height = height, dpi = 180)
  message_time("Wrote ", pdf_path)
  message_time("Wrote ", tif_path)
  message_time("Wrote ", png_path)
  c(pdf = pdf_path, tif = tif_path, png = png_path)
}

draw_heatmap_object <- function(x) {
  if (requireNamespace("ComplexHeatmap", quietly = TRUE) && inherits(x, "Heatmap")) {
    ComplexHeatmap::draw(x, merge_legends = TRUE)
  } else if (requireNamespace("ComplexHeatmap", quietly = TRUE) && inherits(x, "HeatmapList")) {
    ComplexHeatmap::draw(x, merge_legends = TRUE)
  } else {
    print(x)
  }
}

save_heatmap <- function(x, basename, width = 10, height = 8, dpi = 600) {
  plot_obj <- if (is.list(x) && "plot" %in% names(x)) x$plot else x
  pdf_path <- file.path(figure_dir, paste0(basename, ".pdf"))
  tif_path <- file.path(figure_dir, paste0(basename, ".tif"))
  png_path <- file.path(figure_dir, paste0(basename, ".png"))

  if (inherits(plot_obj, c("ggplot", "patchwork"))) {
    ggsave(pdf_path, plot = plot_obj, width = width, height = height, device = cairo_pdf)
    ggsave(tif_path, plot = plot_obj, width = width, height = height, dpi = dpi, device = "tiff", compression = "lzw")
    ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = 180)
  } else {
    grDevices::cairo_pdf(pdf_path, width = width, height = height)
    draw_heatmap_object(plot_obj)
    grDevices::dev.off()
    grDevices::tiff(tif_path, width = width, height = height, units = "in", res = dpi, compression = "lzw")
    draw_heatmap_object(plot_obj)
    grDevices::dev.off()
    grDevices::png(png_path, width = width, height = height, units = "in", res = 180)
    draw_heatmap_object(plot_obj)
    grDevices::dev.off()
  }
  message_time("Wrote ", pdf_path)
  message_time("Wrote ", tif_path)
  message_time("Wrote ", png_path)
  c(pdf = pdf_path, tif = tif_path, png = png_path)
}

ensure_processed <- function(srt, nfeatures = 3000, npcs = 30, prefix = "proc") {
  DefaultAssay(srt) <- "RNA"
  srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  nfeatures <- min(nfeatures, nrow(srt))
  srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  vars <- VariableFeatures(srt)
  if (length(vars) == 0) {
    vars <- rownames(srt)
  }
  srt <- ScaleData(srt, features = vars, verbose = FALSE)
  npcs <- min(npcs, ncol(srt) - 1, length(vars) - 1)
  if (npcs >= 2) {
    srt <- RunPCA(srt, features = vars, npcs = npcs, reduction.name = paste0(prefix, "_pca"), verbose = FALSE)
  }
  srt
}

sanitize_label <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

available_features <- function(features, srt, max_n = Inf) {
  features <- unique(features[!is.na(features) & nzchar(features)])
  features <- features[features %in% rownames(srt)]
  if (is.finite(max_n)) {
    features <- head(features, max_n)
  }
  features
}

marker_features_by_cluster <- function(markers, max_per_cluster = 5, max_total = 80) {
  markers <- markers[!technical_gene(markers$gene), , drop = FALSE]
  markers <- markers[order(markers$cluster, markers$rank), , drop = FALSE]
  by_cluster <- split(markers$gene, markers$cluster)
  unique(head(unlist(lapply(by_cluster, head, max_per_cluster), use.names = FALSE), max_total))
}

marker_features_by_global <- function(markers, annotation_map, max_per_group = 5, max_total = 90) {
  if (!all(c("cluster", "proposed_global_annotation") %in% colnames(annotation_map))) {
    return(character())
  }
  markers <- merge(markers, annotation_map[, c("cluster", "proposed_global_annotation")], by = "cluster")
  markers <- markers[!technical_gene(markers$gene), , drop = FALSE]
  markers <- markers[order(markers$proposed_global_annotation, markers$rank), , drop = FALSE]
  by_group <- split(markers$gene, markers$proposed_global_annotation)
  unique(head(unlist(lapply(by_group, head, max_per_group), use.names = FALSE), max_total))
}

write_matrix <- function(mat, path, row_name = "row") {
  out <- as.data.frame(mat, check.names = FALSE)
  out <- cbind(setNames(data.frame(rownames(out), check.names = FALSE), row_name), out)
  write_tsv(out, path)
}

save_cellcor_outputs <- function(result, basename, width = 10, height = 8) {
  save_heatmap(result, basename, width = width, height = height)
  if (is.list(result) && "simil_matrix" %in% names(result)) {
    write_matrix(
      result$simil_matrix,
      file.path(table_dir, paste0(basename, "_similarity_matrix.tsv")),
      row_name = "query_group"
    )
  }
}

lineage_subsets <- list(
  t_nk = list(
    label = "T/NK",
    pattern = "^T/NK",
    start = "T/NK | ref:CD4_C1",
    end = NULL,
    markers = c("IL7R", "LTB", "CCR7", "CD3D", "CD8A", "CD8B", "GZMK", "GZMB", "NKG7", "GNLY", "PRF1", "KLRD1")
  ),
  myeloid_tam = list(
    label = "Myeloid/TAM",
    pattern = "^Myeloid/TAM",
    start = "Myeloid/TAM | ref:Classical Mono_C1",
    end = NULL,
    markers = c("LYZ", "S100A8", "S100A9", "LST1", "FCN1", "CST3", "AIF1", "CD68", "C1QA", "C1QB", "C1QC", "MRC1")
  ),
  b_plasma = list(
    label = "B/Plasma",
    pattern = "^B/Plasma",
    start = "B/Plasma | ref:B cell_C3",
    end = "B/Plasma | ref:Plasma_C0",
    markers = c("MS4A1", "CD79A", "CD79B", "CD74", "MZB1", "JCHAIN", "IGHG1", "IRF7", "GZMB")
  ),
  epithelial_tumor = list(
    label = "Epithelial/tumor-related",
    pattern = "^Epithelial|^Tumor|^Malignant",
    start = NULL,
    end = NULL,
    markers = c("EPCAM", "KRT8", "KRT18", "KRT19", "MUC1", "MKI67", "TOP2A")
  ),
  endothelial_fibroblast = list(
    label = "Endothelial/Fibroblast",
    pattern = "^Endothelial/Fibroblast",
    start = NULL,
    end = NULL,
    markers = c("PECAM1", "VWF", "COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "RGS5")
  )
)

run_trajectory_subset <- function(srt, subset_id, spec) {
  cells <- rownames(srt@meta.data)[grepl(spec$pattern, as.character(srt$global_annotation))]
  if (length(cells) == 0) {
    record_status(paste0("trajectory_", subset_id), "skipped", "No cells matched the lineage pattern")
    return(NULL)
  }
  sub <- subset(srt, cells = cells)
  sub$global_annotation <- droplevels(factor(sub$global_annotation))
  group_counts <- table(sub$global_annotation)
  valid_groups <- names(group_counts[group_counts >= 30])
  if (length(valid_groups) < 3) {
    record_status(
      paste0("trajectory_", subset_id),
      "skipped",
      paste0("Fewer than 3 groups with at least 30 cells: ", paste(names(group_counts), group_counts, sep = "=", collapse = "; "))
    )
    return(NULL)
  }
  sub <- subset(sub, cells = rownames(sub@meta.data)[sub$global_annotation %in% valid_groups])
  sub$global_annotation <- droplevels(factor(sub$global_annotation))
  message_time("Running trajectory subset ", subset_id, " with ", ncol(sub), " cells and ", length(valid_groups), " groups")

  sub <- ensure_processed(sub, nfeatures = min(2500, nrow(sub)), npcs = 30, prefix = subset_id)
  pca_key <- paste0(subset_id, "_pca")
  dims_use <- seq_len(min(20, ncol(Embeddings(sub, pca_key))))
  sub <- FindNeighbors(sub, reduction = pca_key, dims = dims_use, k.param = 20, verbose = FALSE)
  sub <- FindClusters(sub, resolution = 0.4, algorithm = 1, verbose = FALSE)
  sub <- RunUMAP(
    sub,
    reduction = pca_key,
    dims = dims_use,
    reduction.name = paste0(subset_id, "_umap"),
    reduction.key = paste0(toupper(subset_id), "UMAP_"),
    seed.use = 20260610,
    verbose = FALSE
  )

  start_group <- spec$start
  if (!is.null(start_group) && !start_group %in% levels(sub$global_annotation)) {
    start_group <- NULL
  }
  end_group <- spec$end
  if (!is.null(end_group) && !end_group %in% levels(sub$global_annotation)) {
    end_group <- NULL
  }

  sub <- RunSlingshot(
    sub,
    group.by = "global_annotation",
    reduction = paste0(subset_id, "_umap"),
    prefix = subset_id,
    start = start_group,
    end = end_group,
    show_plot = FALSE,
    seed = 20260610,
    verbose = FALSE
  )

  lineage_cols <- grep(paste0("^", subset_id, "_Lineage"), colnames(sub@meta.data), value = TRUE)
  if (length(lineage_cols) == 0) {
    stop("Slingshot did not produce pseudotime columns for ", subset_id)
  }
  lineage_summary <- do.call(rbind, lapply(lineage_cols, function(col) {
    data.frame(
      subset = subset_id,
      lineage = col,
      cells_with_pseudotime = sum(!is.na(sub@meta.data[[col]])),
      min_pseudotime = suppressWarnings(min(sub@meta.data[[col]], na.rm = TRUE)),
      max_pseudotime = suppressWarnings(max(sub@meta.data[[col]], na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  write_tsv(lineage_summary, file.path(table_dir, paste0("GSE234129_scop_extended_", subset_id, "_pseudotime_summary.tsv")))

  saveRDS(sub, file.path(object_dir, paste0("GSE234129_scop_extended_", subset_id, "_slingshot_seurat.rds")))

  p_lineage <- CellDimPlot(
    sub,
    group.by = "global_annotation",
    reduction = paste0(subset_id, "_umap"),
    lineages = lineage_cols,
    lineages_span = 0.1,
    pt.size = 0.18,
    label = TRUE,
    label_insitu = TRUE,
    theme_use = "theme_blank",
    title = paste0(spec$label, " Slingshot lineages")
  )
  save_plot(p_lineage, paste0("GSE234129_scop_extended_", subset_id, "_slingshot_lineages"), width = 9, height = 7)

  feature_cols <- lineage_cols[seq_len(min(length(lineage_cols), 4))]
  p_pt <- FeatureDimPlot(
    sub,
    features = feature_cols,
    reduction = paste0(subset_id, "_umap"),
    pt.size = 0.18,
    theme_use = "theme_blank",
    title = paste0(spec$label, " pseudotime"),
    ncol = 2
  )
  save_plot(p_pt, paste0("GSE234129_scop_extended_", subset_id, "_pseudotime_umap"), width = 10, height = 8)

  for (time_key in feature_cols) {
    p_proj <- PseudotimeProjectionPlot(
      sub,
      reduction = paste0(subset_id, "_umap"),
      time_key = time_key,
      group.by = "global_annotation",
      plot_type = "stream",
      method = "knn",
      k = 30,
      density = 1.2,
      pt.size = 0.4,
      pt.alpha = 0.25,
      theme_use = "theme_blank",
      title = paste0(spec$label, " ", time_key, " projection")
    )
    save_plot(p_proj, paste0("GSE234129_scop_extended_", subset_id, "_", time_key, "_projection"), width = 9, height = 7)
  }

  dyn_features <- available_features(c(spec$markers, VariableFeatures(sub)), sub, max_n = 60)
  if (length(dyn_features) >= 5) {
    tryCatch(
      {
        dyn <- DynamicHeatmap(
          sub,
          lineages = lineage_cols[1],
          features = dyn_features,
          layer = "data",
          exp_method = "zscore",
          show_row_names = TRUE,
          cluster_rows = TRUE,
          nlabel = 0,
          cores = 1,
          seed = 20260610,
          height = 7,
          width = 10
        )
        save_heatmap(dyn, paste0("GSE234129_scop_extended_", subset_id, "_dynamic_heatmap"), width = 11, height = 8)
        record_status(paste0("trajectory_", subset_id, "_dynamic_heatmap"), "completed", "")
      },
      error = function(e) {
        record_status(paste0("trajectory_", subset_id, "_dynamic_heatmap"), "skipped_or_failed", conditionMessage(e))
        message_time("Skipped/failed dynamic heatmap for ", subset_id, ": ", conditionMessage(e))
      }
    )
  } else {
    record_status(paste0("trajectory_", subset_id, "_dynamic_heatmap"), "skipped", "Fewer than 5 available dynamic features")
  }

  if (any(lineage_summary$cells_with_pseudotime >= 100)) {
    record_status(paste0("trajectory_", subset_id), "completed", paste(lineage_cols, collapse = ", "))
  } else {
    record_status(paste0("trajectory_", subset_id), "skipped_or_failed", "No lineage covered at least 100 cells")
  }
  sub
}

message_time("Reading QC and raw plotting objects")
if (!file.exists(qc_object_path)) {
  stop("Missing QC scop object: ", qc_object_path)
}
if (!file.exists(raw_object_path)) {
  stop("Missing raw scop object: ", raw_object_path)
}
srt <- readRDS(qc_object_path)
raw_srt <- readRDS(raw_object_path)
DefaultAssay(srt) <- "RNA"
DefaultAssay(raw_srt) <- "RNA"

required_meta <- c("patient", "sample", "celltype", "global_annotation", "leiden_0.5")
missing_meta <- setdiff(required_meta, colnames(srt@meta.data))
if (length(missing_meta) > 0) {
  stop("QC object is missing metadata columns: ", paste(missing_meta, collapse = ", "))
}
if (!"umap" %in% names(srt@reductions)) {
  stop("QC object is missing the umap reduction")
}

message_time("Normalizing QC and raw objects for heatmaps and correlations")
srt <- ensure_processed(srt, nfeatures = 3000, npcs = 30, prefix = "qc")
raw_srt <- ensure_processed(raw_srt, nfeatures = 3000, npcs = 30, prefix = "raw")

markers <- read_tsv(marker_path)
cluster_annotation <- read_tsv(cluster_annotation_path)
markers$cluster <- as.character(markers$cluster)
cluster_annotation$cluster <- as.character(cluster_annotation$cluster)

present_classic_markers <- available_features(classic_markers, srt, max_n = 60)
cluster_marker_features <- available_features(marker_features_by_cluster(markers, max_per_cluster = 4, max_total = 80), srt)
global_marker_features <- available_features(marker_features_by_global(markers, cluster_annotation, max_per_group = 5, max_total = 90), srt)
cor_features <- available_features(unique(c(present_classic_markers, global_marker_features, cluster_marker_features)), srt, max_n = 120)

feature_set_table <- data.frame(
  feature_set = c(
    rep("classic_markers", length(present_classic_markers)),
    rep("leiden_top_markers", length(cluster_marker_features)),
    rep("global_annotation_top_markers", length(global_marker_features)),
    rep("correlation_features", length(cor_features))
  ),
  gene = c(present_classic_markers, cluster_marker_features, global_marker_features, cor_features),
  stringsAsFactors = FALSE
)
write_tsv(feature_set_table, file.path(table_dir, "GSE234129_scop_extended_feature_sets.tsv"))

run_step("group_heatmap_classic_markers", {
  ht <- GroupHeatmap(
    srt,
    features = present_classic_markers,
    group.by = c("global_annotation", "leiden_0.5"),
    layer = "data",
    exp_method = "zscore",
    heatmap_palette = "Spectral",
    show_row_names = TRUE,
    row_names_side = "left",
    add_dot = TRUE,
    add_reticle = TRUE,
    nlabel = 0,
    width = 10,
    height = 8
  )
  save_heatmap(ht, "GSE234129_scop_extended_group_heatmap_classic_markers", width = 12, height = 9)
})

run_step("group_heatmap_global_top_markers", {
  ht <- GroupHeatmap(
    srt,
    features = global_marker_features,
    group.by = "global_annotation",
    grouping.var = "sample",
    layer = "data",
    exp_method = "zscore",
    heatmap_palette = "Spectral",
    show_row_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    add_dot = TRUE,
    nlabel = 0,
    width = 11,
    height = 9
  )
  save_heatmap(ht, "GSE234129_scop_extended_group_heatmap_global_annotation_top_markers", width = 13, height = 10)
})

run_step("group_heatmap_leiden_top_markers", {
  ht <- GroupHeatmap(
    srt,
    features = cluster_marker_features,
    group.by = "leiden_0.5",
    grouping.var = "sample",
    layer = "data",
    exp_method = "zscore",
    heatmap_palette = "Spectral",
    show_row_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    add_dot = TRUE,
    nlabel = 0,
    width = 11,
    height = 9
  )
  save_heatmap(ht, "GSE234129_scop_extended_group_heatmap_leiden05_top_markers", width = 13, height = 10)
})

run_step("feature_heatmap_classic_markers_sampled", {
  set.seed(20260610)
  cells_use <- unlist(tapply(colnames(srt), srt$global_annotation, function(x) sample(x, min(length(x), 120))))
  ht <- FeatureHeatmap(
    srt,
    features = present_classic_markers,
    cells = cells_use,
    group.by = "global_annotation",
    layer = "data",
    exp_method = "zscore",
    show_row_names = TRUE,
    show_column_names = FALSE,
    heatmap_palette = "Spectral",
    nlabel = 0,
    width = 11,
    height = 8
  )
  save_heatmap(ht, "GSE234129_scop_extended_feature_heatmap_classic_markers_sampled", width = 12, height = 9)
})

run_step("cellcor_leiden_vs_global_annotation", {
  res <- CellCorHeatmap(
    srt_query = srt,
    query_group = "leiden_0.5",
    ref_group = "global_annotation",
    features = cor_features,
    query_collapsing = TRUE,
    ref_collapsing = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    heatmap_palette = "RdBu",
    nlabel = 0,
    seed = 20260610,
    verbose = FALSE
  )
  save_cellcor_outputs(res, "GSE234129_scop_extended_cellcor_leiden05_vs_global_annotation", width = 12, height = 10)
})

run_step("cellcor_celltype_vs_global_annotation", {
  res <- CellCorHeatmap(
    srt_query = srt,
    query_group = "celltype",
    ref_group = "global_annotation",
    features = cor_features,
    query_collapsing = TRUE,
    ref_collapsing = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    heatmap_palette = "RdBu",
    nlabel = 0,
    seed = 20260610,
    verbose = FALSE
  )
  save_cellcor_outputs(res, "GSE234129_scop_extended_cellcor_celltype_vs_global_annotation", width = 13, height = 12)
})

run_step("cellcor_sample_similarity", {
  res <- CellCorHeatmap(
    srt_query = srt,
    query_group = "sample",
    ref_group = "sample",
    features = cor_features,
    query_collapsing = TRUE,
    ref_collapsing = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    heatmap_palette = "RdBu",
    nlabel = 0,
    seed = 20260610,
    verbose = FALSE
  )
  save_cellcor_outputs(res, "GSE234129_scop_extended_cellcor_sample_similarity", width = 10, height = 9)
})

run_step("cellcor_raw_celltype_vs_qc_global_annotation", {
  raw_features <- available_features(cor_features, raw_srt, max_n = 120)
  res <- CellCorHeatmap(
    srt_query = srt,
    srt_ref = raw_srt,
    query_group = "global_annotation",
    ref_group = "celltype",
    features = raw_features,
    query_collapsing = TRUE,
    ref_collapsing = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    heatmap_palette = "RdBu",
    nlabel = 0,
    seed = 20260610,
    verbose = FALSE
  )
  save_cellcor_outputs(res, "GSE234129_scop_extended_cellcor_qc_global_annotation_vs_raw_celltype", width = 13, height = 10)
})

run_step("cellcor_raw_leiden_vs_qc_leiden", {
  raw_features <- available_features(cor_features, raw_srt, max_n = 120)
  res <- CellCorHeatmap(
    srt_query = srt,
    srt_ref = raw_srt,
    query_group = "leiden_0.5",
    ref_group = "leiden_0.5",
    features = raw_features,
    query_collapsing = TRUE,
    ref_collapsing = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    heatmap_palette = "RdBu",
    nlabel = 0,
    seed = 20260610,
    verbose = FALSE
  )
  save_cellcor_outputs(res, "GSE234129_scop_extended_cellcor_qc_leiden05_vs_raw_leiden05", width = 11, height = 10)
})

trajectory_objects <- list()
for (subset_id in names(lineage_subsets)) {
  run_step(paste0("trajectory_wrapper_", subset_id), {
    trajectory_objects[[subset_id]] <- run_trajectory_subset(srt, subset_id, lineage_subsets[[subset_id]])
  })
}

write_tsv(step_results, summary_path)

figure_files <- list.files(figure_dir, pattern = "\\.(pdf|tif)$", full.names = TRUE)
table_files <- list.files(table_dir, pattern = "\\.tsv$", full.names = TRUE)
object_files <- list.files(object_dir, pattern = "\\.rds$", full.names = TRUE)
optional_pkgs <- c("monocle3", "tradeSeq", "palantir")
optional_status <- vapply(optional_pkgs, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else "NOT_INSTALLED"
}, character(1))

report_lines <- c(
  "# GSE234129 scop Extended Analysis Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  paste0("- QC scop object: `", qc_object_path, "`"),
  paste0("- Raw scop object: `", raw_object_path, "`"),
  paste0("- Marker table: `", marker_path, "`"),
  paste0("- Cluster annotation table: `", cluster_annotation_path, "`"),
  "",
  "## Object Summary",
  "",
  paste0("- QC dimensions: ", nrow(srt), " genes x ", ncol(srt), " cells"),
  paste0("- Raw dimensions: ", nrow(raw_srt), " genes x ", ncol(raw_srt), " cells"),
  paste0("- QC global annotations: ", length(unique(srt$global_annotation))),
  paste0("- QC Leiden 0.5 clusters: ", length(unique(srt$leiden_0.5))),
  paste0("- Raw Leiden 0.5 clusters: ", length(unique(raw_srt$leiden_0.5))),
  "",
  "## Completed Modules",
  "",
  "- GroupHeatmap for classic markers, global annotation top markers, and Leiden top markers.",
  "- FeatureHeatmap for sampled classic marker expression.",
  "- CellCorHeatmap within the QC object and between QC/raw references.",
  "- Slingshot trajectory inference for eligible lineage subsets.",
  "- DynamicHeatmap attempted for each successful lineage with available marker/HVF features.",
  "",
  "## Trajectory Dependency Status",
  "",
  paste0("- slingshot: `", as.character(utils::packageVersion("slingshot")), "`"),
  paste0("- monocle3: `", optional_status[["monocle3"]], "`"),
  paste0("- tradeSeq: `", optional_status[["tradeSeq"]], "`"),
  paste0("- palantir: `", optional_status[["palantir"]], "`"),
  "",
  "## Outputs",
  "",
  paste0("- Figure directory: `", figure_dir, "`"),
  paste0("- Table directory: `", table_dir, "`"),
  paste0("- Local trajectory object directory: `", object_dir, "`"),
  paste0("- Step summary: `", summary_path, "`"),
  paste0("- Figure files written (PDF/TIF): ", length(figure_files)),
  paste0("- Table files written: ", length(table_files)),
  paste0("- Local trajectory RDS objects written: ", length(object_files)),
  "",
  "## Step Status",
  "",
  paste(capture.output(print(step_results, row.names = FALSE)), collapse = "\n"),
  "",
  "## Session Info",
  "",
  "```text",
  capture.output(sessionInfo()),
  "```"
)
writeLines(report_lines, report_path)
message_time("Wrote ", report_path)
message_time("Extended scop analysis complete")
