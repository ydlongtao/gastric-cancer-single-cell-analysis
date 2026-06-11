suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(SingleCellExperiment)
  library(ggplot2)
})

if (!requireNamespace("scrapper", quietly = TRUE)) {
  stop(
    "Missing required Bioconductor package 'scrapper'. Install with: ",
    "BiocManager::install('scrapper', ask = FALSE, update = FALSE)"
  )
}

set.seed(20260611)

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

original_seurat_path <- file.path(project_dir, "Datasets/GSE234129/processed/GSE234129_seurat.rds")
qc_metadata_path <- file.path(project_dir, "results/GSE234129/scdesign3/tables/GSE234129_scdesign3_input_metadata.tsv")
simulated_seurat_path <- file.path(project_dir, "results/GSE234129/scdesign3/objects/GSE234129_scdesign3_simulated_4x_200features_hvg_markers_seurat.rds")

out_object_path <- file.path(project_dir, "results/GSE234129/objects/GSE234129_raw_simulated_200gene_blind_annotated_seurat.rds")
table_dir <- file.path(project_dir, "results/GSE234129/tables/blind_annotation")
figure_dir <- file.path(project_dir, "results/GSE234129/figures/blind_annotation")
report_path <- file.path(project_dir, "results/GSE234129/reports/GSE234129_raw_simulated_200gene_blind_annotation_comparison_report.md")
dir.create(dirname(out_object_path), recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

cluster_key <- "blind_leiden_0.6"

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

marker_panels <- list(
  "Epithelial/Tumor" = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1", "CLDN4", "TACSTD2"),
  "T cell" = c("CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CD4", "CD8A", "CD8B"),
  "NK cell" = c("NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "KLRB1", "XCL1", "XCL2"),
  "B cell" = c("MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA"),
  "Plasma cell" = c("MZB1", "JCHAIN", "IGHG1", "IGHG3", "IGKC", "XBP1", "SDC1"),
  "Myeloid/TAM" = c("LYZ", "LST1", "CST3", "AIF1", "FCER1G", "CD68", "MRC1", "C1QA", "C1QB", "SPP1"),
  "Mast cell" = c("TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"),
  "Endothelial" = c("PECAM1", "VWF", "KDR", "CLDN5", "ENG", "ESAM"),
  "Fibroblast/Pericyte" = c("COL1A1", "COL1A2", "DCN", "LUM", "ACTA2", "RGS5", "MCAM", "PDGFRB")
)

hpca_to_lineage <- list(
  "T cell" = c("T_cells", "T_cell", "T.cells", "T cells", "T_cell:CD4+", "T_cell:CD8+"),
  "NK cell" = c("NK_cell", "NK.cells", "NK cells"),
  "B cell" = c("B_cell", "B.cells", "B cells"),
  "Plasma cell" = c("Plasma_cell", "Plasma.cells", "Plasma cells"),
  "Myeloid/TAM" = c("Monocyte", "Macrophage", "Dendritic_cell", "DC", "Neutrophils"),
  "Mast cell" = c("Mast_cell", "Mast.cells", "Mast cells"),
  "Endothelial" = c("Endothelial_cells", "Endothelial.cell", "Endothelial cells"),
  "Fibroblast/Pericyte" = c("Fibroblasts", "Smooth_muscle_cells", "MSC", "Pericytes"),
  "Epithelial/Tumor" = c("Epithelial_cells", "Epithelial.cell", "Epithelial cells")
)

message_time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

save_plot <- function(plot, basename, width = 8, height = 6, dpi = 600) {
  pdf_path <- file.path(figure_dir, paste0(basename, ".pdf"))
  tif_path <- file.path(figure_dir, paste0(basename, ".tif"))
  ggsave(pdf_path, plot = plot, width = width, height = height, device = cairo_pdf)
  ggsave(tif_path, plot = plot, width = width, height = height, dpi = dpi, device = "tiff", compression = "lzw")
  message_time("Wrote ", pdf_path)
  message_time("Wrote ", tif_path)
}

read_required <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " not found: ", path)
  }
  invisible(path)
}

normalize_holdout <- function(x) {
  out <- as.character(x)
  out[is.na(out) | out == ""] <- "Unassigned"
  out
}

broad_from_previous_label <- function(x) {
  x <- as.character(x)
  out <- rep("Other/Unknown", length(x))
  out[grepl("CD4|CD8|T/NK|NKT|DNT|NK", x, ignore.case = TRUE)] <- "T/NK"
  out[grepl("B cell|Plasma|pDC", x, ignore.case = TRUE)] <- "B/Plasma"
  out[grepl("TAM|Mono|Macrophage|Myeloid|Mast", x, ignore.case = TRUE)] <- "Myeloid/TAM"
  out[grepl("Endothelial|VSMC|Fibro|Pericyte|COL1|DCN", x, ignore.case = TRUE)] <- "Stromal"
  out[grepl("Epithelial|Tumor|Cancer|KRT|EPCAM", x, ignore.case = TRUE)] <- "Epithelial/Tumor"
  out
}

blind_to_broad <- function(x) {
  out <- as.character(x)
  out[out %in% c("T cell", "NK cell")] <- "T/NK"
  out[out %in% c("B cell", "Plasma cell")] <- "B/Plasma"
  out[out %in% c("Endothelial", "Fibroblast/Pericyte")] <- "Stromal"
  out
}

mode_fraction <- function(values) {
  tab <- sort(table(values), decreasing = TRUE)
  if (length(tab) == 0) {
    return(list(label = NA_character_, fraction = NA_real_))
  }
  list(label = names(tab)[[1]], fraction = as.numeric(tab[[1]]) / sum(tab))
}

contingency_long <- function(left, right, left_name, right_name) {
  tab <- as.data.frame.matrix(table(left, right))
  tab[[left_name]] <- rownames(tab)
  out <- reshape(
    tab,
    varying = setdiff(colnames(tab), left_name),
    v.names = "n_cells",
    timevar = right_name,
    times = setdiff(colnames(tab), left_name),
    direction = "long"
  )
  rownames(out) <- NULL
  out <- out[out$n_cells > 0, c(left_name, right_name, "n_cells")]
  totals <- aggregate(out$n_cells, by = list(out[[left_name]]), FUN = sum)
  colnames(totals) <- c(left_name, "blind_label_total")
  out <- merge(out, totals, by = left_name, all.x = TRUE)
  out$fraction_within_blind_label <- out$n_cells / out$blind_label_total
  out[order(out[[left_name]], -out$n_cells), , drop = FALSE]
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n < 2) {
    return(NA_real_)
  }
  comb2 <- function(v) v * (v - 1) / 2
  sum_nij <- sum(comb2(tab))
  sum_ai <- sum(comb2(rowSums(tab)))
  sum_bj <- sum(comb2(colSums(tab)))
  total <- comb2(n)
  expected <- sum_ai * sum_bj / total
  max_index <- 0.5 * (sum_ai + sum_bj)
  denom <- max_index - expected
  if (denom == 0) {
    return(NA_real_)
  }
  (sum_nij - expected) / denom
}

normalized_mutual_information <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  if (n == 0) {
    return(NA_real_)
  }
  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- pxy > 0
  mi <- sum(pxy[nz] * log(pxy[nz] / outer(px, py)[nz]))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx + hy == 0) {
    return(NA_real_)
  }
  2 * mi / (hx + hy)
}

purity_summary <- function(pred, ref, metric_name) {
  split_ref <- split(ref, pred)
  rows <- lapply(names(split_ref), function(label) {
    mf <- mode_fraction(split_ref[[label]])
    data.frame(
      comparison = metric_name,
      blind_annotation = label,
      dominant_holdout_label = mf$label,
      dominant_holdout_fraction = mf$fraction,
      n_cells = length(split_ref[[label]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

heatmap_plot <- function(table_df, x_col, y_col, title) {
  plot_df <- table_df
  plot_df[[x_col]] <- factor(plot_df[[x_col]], levels = unique(plot_df[[x_col]]))
  plot_df[[y_col]] <- factor(plot_df[[y_col]], levels = rev(unique(plot_df[[y_col]])))
  ggplot(plot_df, aes(x = .data[[x_col]], y = .data[[y_col]], fill = fraction_within_blind_label)) +
    geom_tile(color = "white", linewidth = 0.15) +
    scale_fill_gradient(low = "#f7fbff", high = "#145da0", limits = c(0, 1), name = "Fraction") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size = 6)) +
    labs(x = x_col, y = y_col, title = title)
}

message_time("Checking inputs")
read_required(original_seurat_path, "Original Seurat object")
read_required(qc_metadata_path, "QC metadata")
read_required(simulated_seurat_path, "Simulated Seurat object")

message_time("Loading objects")
original <- readRDS(original_seurat_path)
simulated <- readRDS(simulated_seurat_path)
qc_metadata <- read.delim(qc_metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!"cell_barcode" %in% colnames(qc_metadata)) {
  stop("QC metadata is missing cell_barcode column: ", qc_metadata_path)
}

qc_cells <- intersect(qc_metadata$cell_barcode, Cells(original))
if (length(qc_cells) == 0) {
  stop("No QC metadata cell barcodes matched the original Seurat object")
}

shared_genes <- intersect(rownames(simulated), rownames(original))
if (length(shared_genes) < 20) {
  stop("Too few shared genes between original and simulated objects: ", length(shared_genes))
}
message_time("Shared genes: ", length(shared_genes), "; QC original cells: ", length(qc_cells), "; simulated cells: ", length(Cells(simulated)))

original_counts <- GetAssayData(original, assay = "RNA", layer = "counts")[shared_genes, qc_cells, drop = FALSE]
simulated_counts <- GetAssayData(simulated, assay = "RNA", layer = "counts")[shared_genes, Cells(simulated), drop = FALSE]
colnames(original_counts) <- paste0("original_", colnames(original_counts))
colnames(simulated_counts) <- paste0("simulated_", colnames(simulated_counts))
combined_counts <- cbind(original_counts, simulated_counts)

original_meta <- original@meta.data[qc_cells, , drop = FALSE]
original_qc_meta <- qc_metadata[match(qc_cells, qc_metadata$cell_barcode), , drop = FALSE]
simulated_meta <- simulated@meta.data[Cells(simulated), , drop = FALSE]

analysis_meta <- rbind(
  data.frame(
    cell_barcode = colnames(original_counts),
    source = "original_qc",
    sample = normalize_holdout(original_meta$sample),
    stringsAsFactors = FALSE
  ),
  data.frame(
    cell_barcode = colnames(simulated_counts),
    source = "scdesign3_simulated",
    sample = normalize_holdout(simulated_meta$sample),
    stringsAsFactors = FALSE
  )
)
rownames(analysis_meta) <- analysis_meta$cell_barcode

holdout_meta <- rbind(
  data.frame(
    cell_barcode = colnames(original_counts),
    holdout_original_celltype = normalize_holdout(original_meta$celltype),
    holdout_previous_global_annotation = normalize_holdout(original_qc_meta$global_annotation),
    holdout_simulated_reference_annotation = NA_character_,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cell_barcode = colnames(simulated_counts),
    holdout_original_celltype = NA_character_,
    holdout_previous_global_annotation = NA_character_,
    holdout_simulated_reference_annotation = normalize_holdout(simulated_meta$global_annotation),
    stringsAsFactors = FALSE
  )
)
rownames(holdout_meta) <- holdout_meta$cell_barcode

message_time("Creating merged Seurat object without holdout labels in analysis metadata")
combined <- CreateSeuratObject(
  counts = combined_counts,
  meta.data = analysis_meta,
  assay = "RNA",
  project = "GSE234129_raw_simulated_200gene_blind"
)
combined$source <- factor(combined$source, levels = c("original_qc", "scdesign3_simulated"))

message_time("Normalizing, clustering, and embedding merged 200-gene object")
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = min(2000, nrow(combined)), verbose = FALSE)
combined <- ScaleData(combined, features = rownames(combined), verbose = FALSE)
combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = min(30, length(VariableFeatures(combined)) - 1), verbose = FALSE)
combined <- FindNeighbors(combined, dims = 1:20, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.6, algorithm = 1, random.seed = 20260611, verbose = FALSE)
combined[[cluster_key]] <- combined$seurat_clusters
combined <- RunUMAP(combined, dims = 1:20, seed.use = 20260611, verbose = FALSE)

message_time("Calculating cluster markers")
Idents(combined) <- cluster_key
cluster_markers <- FindAllMarkers(
  combined,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.10,
  test.use = "wilcox",
  verbose = FALSE
)
write.table(
  cluster_markers,
  file.path(table_dir, "GSE234129_raw_simulated_200gene_cluster_markers.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message_time("Scoring searched marker panels")
present_panel_genes <- unique(unlist(marker_panels))
present_panel_genes <- present_panel_genes[present_panel_genes %in% rownames(combined)]
if (length(present_panel_genes) == 0) {
  stop("No marker panel genes are present in the merged 200-gene object")
}

log_data <- GetAssayData(combined, assay = "RNA", layer = "data")
cluster_ids <- as.character(combined[[cluster_key]][, 1])
marker_score_rows <- list()
for (panel_name in names(marker_panels)) {
  genes <- intersect(marker_panels[[panel_name]], rownames(log_data))
  if (length(genes) == 0) {
    next
  }
  values <- Matrix::colMeans(log_data[genes, , drop = FALSE])
  score_df <- data.frame(cluster = cluster_ids, score = as.numeric(values))
  cluster_scores <- aggregate(score ~ cluster, score_df, mean)
  cluster_scores$marker_panel <- panel_name
  cluster_scores$genes_present <- paste(genes, collapse = ",")
  marker_score_rows[[panel_name]] <- cluster_scores[, c("cluster", "marker_panel", "score", "genes_present")]
}
marker_scores <- do.call(rbind, marker_score_rows)
rownames(marker_scores) <- NULL
write.table(
  marker_scores,
  file.path(table_dir, "GSE234129_raw_simulated_200gene_marker_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message_time("Running cluster-level SingleR with celldex HPCA where shared genes permit")
singler_table <- data.frame(
  cluster = sort(unique(cluster_ids)),
  singler_label = NA_character_,
  singler_pruned_label = NA_character_,
  singler_lineage = NA_character_,
  singler_status = "not_run",
  stringsAsFactors = FALSE
)
singler_error <- NA_character_
tryCatch({
  hpca <- celldex::HumanPrimaryCellAtlasData()
  hpca_genes <- intersect(rownames(hpca), rownames(log_data))
  if (length(hpca_genes) < 20) {
    stop("HPCA shares fewer than 20 genes with the 200-gene merged object")
  }
  sce <- SingleCellExperiment(list(logcounts = as.matrix(log_data[hpca_genes, , drop = FALSE])))
  pred <- SingleR(
    test = sce,
    ref = hpca[hpca_genes, ],
    labels = hpca$label.fine,
    clusters = cluster_ids
  )
  singler_table <- data.frame(
    cluster = rownames(pred),
    singler_label = as.character(pred$labels),
    singler_pruned_label = normalize_holdout(pred$pruned.labels),
    singler_lineage = NA_character_,
    singler_status = "ok",
    stringsAsFactors = FALSE
  )
  for (lineage in names(hpca_to_lineage)) {
    hits <- grepl(paste(hpca_to_lineage[[lineage]], collapse = "|"), singler_table$singler_label, ignore.case = TRUE)
    singler_table$singler_lineage[hits & is.na(singler_table$singler_lineage)] <- lineage
  }
  singler_table$singler_lineage[is.na(singler_table$singler_lineage)] <- "Unmapped HPCA label"
}, error = function(e) {
  singler_error <<- conditionMessage(e)
  singler_table$singler_status <<- "failed"
})

write.table(
  singler_table,
  file.path(table_dir, "GSE234129_raw_simulated_200gene_singler_hpca_cluster_predictions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message_time("Assigning blind cluster labels by marker and SingleR consensus")
cluster_levels <- sort(unique(cluster_ids))
annotation_rows <- lapply(cluster_levels, function(cl) {
  scores <- marker_scores[marker_scores$cluster == cl, , drop = FALSE]
  scores <- scores[order(-scores$score), , drop = FALSE]
  best_panel <- scores$marker_panel[[1]]
  best_score <- scores$score[[1]]
  second_score <- if (nrow(scores) >= 2) scores$score[[2]] else 0
  top_genes <- cluster_markers[as.character(cluster_markers$cluster) == cl, , drop = FALSE]
  top_genes <- head(top_genes$gene, 10)
  sr <- singler_table[singler_table$cluster == cl, , drop = FALSE]
  sr_lineage <- if (nrow(sr) == 1) sr$singler_lineage[[1]] else NA_character_
  sr_label <- if (nrow(sr) == 1) sr$singler_label[[1]] else NA_character_
  marker_margin <- best_score - second_score
  consensus <- !is.na(sr_lineage) && sr_lineage == best_panel
  label <- best_panel
  confidence <- if (isTRUE(consensus) && marker_margin >= 0.10 && best_score > 0.15) {
    "high"
  } else if ((best_score > 0.12 && marker_margin >= 0.04) || isTRUE(consensus)) {
    "medium"
  } else {
    "low"
  }
  n_cells <- sum(cluster_ids == cl)
  source_mix <- prop.table(table(combined$source[cluster_ids == cl]))
  data.frame(
    cluster = cl,
    n_cells = n_cells,
    blind_annotation = label,
    blind_annotation_confidence = confidence,
    best_marker_score = best_score,
    second_marker_score = second_score,
    marker_score_margin = marker_margin,
    singler_hpca_label = normalize_holdout(sr_label),
    singler_mapped_lineage = normalize_holdout(sr_lineage),
    source_original_fraction = as.numeric(source_mix["original_qc"] %||% 0),
    source_simulated_fraction = as.numeric(source_mix["scdesign3_simulated"] %||% 0),
    blind_annotation_evidence = paste0(
      "markers=", paste(top_genes, collapse = ","),
      "; marker_panel=", best_panel,
      "; SingleR=", normalize_holdout(sr_label),
      "; SingleR_lineage=", normalize_holdout(sr_lineage)
    ),
    stringsAsFactors = FALSE
  )
})
cluster_annotation <- do.call(rbind, annotation_rows)
write.table(
  cluster_annotation,
  file.path(table_dir, "GSE234129_raw_simulated_200gene_cluster_annotation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

label_map <- setNames(cluster_annotation$blind_annotation, cluster_annotation$cluster)
confidence_map <- setNames(cluster_annotation$blind_annotation_confidence, cluster_annotation$cluster)
evidence_map <- setNames(cluster_annotation$blind_annotation_evidence, cluster_annotation$cluster)
combined$blind_annotation <- factor(unname(label_map[cluster_ids]))
combined$blind_annotation_confidence <- factor(unname(confidence_map[cluster_ids]), levels = c("high", "medium", "low"))
combined$blind_annotation_evidence <- unname(evidence_map[cluster_ids])

message_time("Reattaching holdout labels after blind annotation")
holdout_meta <- holdout_meta[Cells(combined), , drop = FALSE]
combined$holdout_original_celltype <- holdout_meta$holdout_original_celltype
combined$holdout_previous_global_annotation <- holdout_meta$holdout_previous_global_annotation
combined$holdout_simulated_reference_annotation <- holdout_meta$holdout_simulated_reference_annotation

meta <- combined@meta.data
original_cells <- meta$source == "original_qc"
simulated_cells <- meta$source == "scdesign3_simulated"

original_comparison <- contingency_long(
  combined$blind_annotation[original_cells],
  combined$holdout_original_celltype[original_cells],
  "blind_annotation",
  "holdout_original_celltype"
)
simulated_comparison <- contingency_long(
  combined$blind_annotation[simulated_cells],
  combined$holdout_simulated_reference_annotation[simulated_cells],
  "blind_annotation",
  "holdout_simulated_reference_annotation"
)
write.table(
  original_comparison,
  file.path(table_dir, "GSE234129_raw_original_blind_vs_holdout_celltype.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  simulated_comparison,
  file.path(table_dir, "GSE234129_simulated_blind_vs_holdout_reference.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

original_broad_holdout <- broad_from_previous_label(combined$holdout_original_celltype[original_cells])
simulated_broad_holdout <- broad_from_previous_label(combined$holdout_simulated_reference_annotation[simulated_cells])
original_broad_blind <- blind_to_broad(combined$blind_annotation[original_cells])
simulated_broad_blind <- blind_to_broad(combined$blind_annotation[simulated_cells])

source_tab <- as.data.frame.matrix(table(combined[[cluster_key]][, 1], combined$source))
source_tab$cluster <- rownames(source_tab)
source_tab$total_cells <- rowSums(source_tab[, setdiff(colnames(source_tab), "cluster"), drop = FALSE])
source_tab$source_balance_score <- apply(source_tab[, c("original_qc", "scdesign3_simulated"), drop = FALSE], 1, function(x) {
  if (sum(x) == 0) {
    return(NA_real_)
  }
  min(x) / sum(x)
})

quality_metrics <- rbind(
  data.frame(metric = "original_blind_vs_celltype_ari", value = adjusted_rand_index(combined$blind_annotation[original_cells], combined$holdout_original_celltype[original_cells])),
  data.frame(metric = "original_blind_vs_celltype_nmi", value = normalized_mutual_information(combined$blind_annotation[original_cells], combined$holdout_original_celltype[original_cells])),
  data.frame(metric = "original_broad_blind_vs_celltype_broad_ari", value = adjusted_rand_index(original_broad_blind, original_broad_holdout)),
  data.frame(metric = "original_broad_blind_vs_celltype_broad_nmi", value = normalized_mutual_information(original_broad_blind, original_broad_holdout)),
  data.frame(metric = "simulated_blind_vs_reference_ari", value = adjusted_rand_index(combined$blind_annotation[simulated_cells], combined$holdout_simulated_reference_annotation[simulated_cells])),
  data.frame(metric = "simulated_blind_vs_reference_nmi", value = normalized_mutual_information(combined$blind_annotation[simulated_cells], combined$holdout_simulated_reference_annotation[simulated_cells])),
  data.frame(metric = "simulated_broad_blind_vs_reference_broad_ari", value = adjusted_rand_index(simulated_broad_blind, simulated_broad_holdout)),
  data.frame(metric = "simulated_broad_blind_vs_reference_broad_nmi", value = normalized_mutual_information(simulated_broad_blind, simulated_broad_holdout)),
  data.frame(metric = "mean_cluster_source_balance_score", value = mean(source_tab$source_balance_score, na.rm = TRUE)),
  data.frame(metric = "median_cluster_source_balance_score", value = median(source_tab$source_balance_score, na.rm = TRUE))
)

purity_table <- rbind(
  purity_summary(combined$blind_annotation[original_cells], combined$holdout_original_celltype[original_cells], "original_celltype"),
  purity_summary(combined$blind_annotation[simulated_cells], combined$holdout_simulated_reference_annotation[simulated_cells], "simulated_reference")
)
write.table(
  quality_metrics,
  file.path(table_dir, "GSE234129_annotation_quality_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  purity_table,
  file.path(table_dir, "GSE234129_raw_simulated_200gene_label_purity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  source_tab,
  file.path(table_dir, "GSE234129_raw_simulated_200gene_cluster_source_mixing.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message_time("Writing figures")
p_source <- DimPlot(combined, reduction = "umap", group.by = "source", pt.size = 0.08) +
  theme_bw() +
  labs(title = "Merged original QC and scDesign3 simulated cells", color = "Source")
save_plot(p_source, "GSE234129_raw_simulated_200gene_umap_source", width = 7, height = 5)

p_blind <- DimPlot(combined, reduction = "umap", group.by = "blind_annotation", pt.size = 0.08, label = TRUE, repel = TRUE) +
  theme_bw() +
  labs(title = "Blind annotation on merged 200-gene object", color = "Blind annotation")
save_plot(p_blind, "GSE234129_raw_simulated_200gene_umap_blind_annotation", width = 8, height = 6)

p_original_heatmap <- heatmap_plot(
  original_comparison,
  "blind_annotation",
  "holdout_original_celltype",
  "Original cells: blind annotation vs held-out original celltype"
)
save_plot(p_original_heatmap, "GSE234129_raw_original_blind_vs_holdout_celltype_heatmap", width = 8.5, height = 7)

p_sim_heatmap <- heatmap_plot(
  simulated_comparison,
  "blind_annotation",
  "holdout_simulated_reference_annotation",
  "Simulated cells: blind annotation vs held-out reference label"
)
save_plot(p_sim_heatmap, "GSE234129_simulated_blind_vs_holdout_reference_heatmap", width = 8.5, height = 7)

dotplot_genes <- unique(unlist(marker_panels))
dotplot_genes <- dotplot_genes[dotplot_genes %in% rownames(combined)]
p_dot <- DotPlot(combined, features = dotplot_genes, group.by = "blind_annotation") +
  theme_bw(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Marker panel expression by blind annotation", x = "Marker gene", y = "Blind annotation")
save_plot(p_dot, "GSE234129_raw_simulated_200gene_marker_dotplot", width = 10, height = 5.5)

message_time("Saving merged object")
combined@misc$blind_annotation_notes <- list(
  date = as.character(Sys.Date()),
  strategy = "Existing 200-gene scDesign3 simulation merged with QC-matched original cells; blind annotation did not use holdout labels.",
  reference_resources = c("SingleR/celldex HPCA", "CellMarker 2.0", "PanglaoDB", "TISCH2", "CellTypist documentation"),
  limitation = "The merged matrix has only 200 genes, limiting fine subtype annotation."
)
saveRDS(combined, out_object_path)

message_time("Writing report")
metrics_md <- paste(capture.output(print(quality_metrics, row.names = FALSE)), collapse = "\n")
cluster_md <- paste(capture.output(print(cluster_annotation[, c("cluster", "n_cells", "blind_annotation", "blind_annotation_confidence", "singler_hpca_label", "source_original_fraction", "source_simulated_fraction")], row.names = FALSE)), collapse = "\n")
singler_note <- if (is.na(singler_error)) {
  "SingleR/celldex HPCA completed at cluster level using genes shared with the 200-gene merged object."
} else {
  paste0("SingleR/celldex HPCA failed and marker-consensus labels were used alone. Error: ", singler_error)
}

report_lines <- c(
  "# GSE234129 Raw + Simulated 200-Gene Blind Annotation Comparison",
  "",
  "## Inputs",
  "",
  paste0("- Original Seurat object: `", original_seurat_path, "`"),
  paste0("- QC metadata: `", qc_metadata_path, "`"),
  paste0("- Simulated Seurat object: `", simulated_seurat_path, "`"),
  "",
  "## Method",
  "",
  paste0("- Merged ", length(qc_cells), " QC-matched original cells and ", ncol(simulated_counts), " scDesign3 simulated cells on ", length(shared_genes), " shared genes."),
  "- The clustering and blind annotation metadata contained only source/sample fields; original `celltype`, previous `global_annotation`, and simulated reference labels were reattached only after `blind_annotation` was assigned.",
  paste0("- ", singler_note),
  "- Marker panels were curated from the searched reference resources: CellMarker 2.0, PanglaoDB, TISCH2, and CellTypist documentation, with gastric cancer/TME-relevant broad lineages.",
  "",
  "## Searched Annotation Resources",
  "",
  "- SingleR/celldex HPCA: https://bioconductor.org/books/release/SingleRBook/introduction.html",
  "- CellMarker 2.0: https://academic.oup.com/nar/article/51/D1/D870/6775381",
  "- PanglaoDB: https://panglaodb.se/",
  "- TISCH2: https://tisch.compbio.cn/documentation/",
  "- CellTypist: https://celltypist.readthedocs.io/",
  "",
  "## Cluster Annotation Summary",
  "",
  "```text",
  cluster_md,
  "```",
  "",
  "## Quantitative Comparison",
  "",
  "```text",
  metrics_md,
  "```",
  "",
  "## Interpretation",
  "",
  "- Agreement should be interpreted primarily at broad lineage level because the selected scDesign3 simulation contains only 200 genes.",
  "- Blind labels that collapse several held-out subtypes are expected when the old labels distinguish fine CD4/CD8/TAM/plasma states unsupported by the restricted gene set.",
  "- Higher-confidence blind labels require both a clear marker-panel margin and agreement with SingleR HPCA where available.",
  "- Low-confidence clusters should not be used for fine biological claims without rerunning scDesign3 with a larger gene set.",
  "",
  "## Outputs",
  "",
  paste0("- Object: `", out_object_path, "`"),
  paste0("- Tables: `", table_dir, "`"),
  paste0("- Figures: `", figure_dir, "`")
)
writeLines(report_lines, report_path)

message_time("Done")
message_time("Object: ", out_object_path)
message_time("Report: ", report_path)
