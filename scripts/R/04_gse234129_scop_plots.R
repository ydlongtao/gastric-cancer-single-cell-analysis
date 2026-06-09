suppressPackageStartupMessages({
  library(Seurat)
  library(scop)
  library(ggplot2)
})

project_dir <- normalizePath(file.path(dirname(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][[1]])), "../.."), mustWork = TRUE)
seurat_path <- file.path(project_dir, "Datasets/GSE234129/processed/GSE234129_seurat.rds")
metadata_path <- file.path(project_dir, "results/GSE234129/tables/GSE234129_scop_plot_metadata.tsv")
out_dir <- file.path(project_dir, "results/GSE234129/figures/scop")
object_out <- file.path(project_dir, "results/GSE234129/objects/GSE234129_scop_plotting_seurat.rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(object_out), recursive = TRUE, showWarnings = FALSE)

save_plot <- function(plot, basename, width = 8, height = 6, dpi = 600) {
  pdf_path <- file.path(out_dir, paste0(basename, ".pdf"))
  tif_path <- file.path(out_dir, paste0(basename, ".tif"))
  ggsave(pdf_path, plot = plot, width = width, height = height, device = cairo_pdf)
  ggsave(
    tif_path,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    device = "tiff",
    compression = "lzw"
  )
  message("Wrote ", pdf_path)
  message("Wrote ", tif_path)
}

srt <- readRDS(seurat_path)
meta <- read.delim(metadata_path, stringsAsFactors = FALSE)
rownames(meta) <- meta$cell_barcode

common_cells <- intersect(colnames(srt), rownames(meta))
if (length(common_cells) == 0) {
  stop("No overlapping cells between Seurat object and exported Scanpy metadata")
}
srt <- subset(srt, cells = common_cells)
meta <- meta[colnames(srt), , drop = FALSE]

srt$patient <- meta$patient
srt$sample <- meta$sample
srt$celltype <- meta$celltype
srt$n_genes_by_counts <- meta$n_genes_by_counts
srt$total_counts <- meta$total_counts
srt$pct_counts_mt <- meta$pct_counts_mt
srt$leiden_0.5 <- factor(meta$leiden_0.5)
srt$global_annotation <- factor(meta$global_annotation)

umap <- as.matrix(meta[, c("UMAP_1", "UMAP_2")])
rownames(umap) <- rownames(meta)
colnames(umap) <- c("UMAP_1", "UMAP_2")
srt[["umap"]] <- CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = DefaultAssay(srt)
)

saveRDS(srt, object_out)
message("Wrote ", object_out)

p_umap_cluster <- CellDimPlot(
  srt,
  group.by = "leiden_0.5",
  reduction = "umap",
  label = TRUE,
  label_insitu = TRUE,
  pt.size = 0.25,
  theme_use = "theme_blank",
  title = "GSE234129 Leiden 0.5 clusters"
)
save_plot(p_umap_cluster, "GSE234129_scop_umap_leiden05", width = 8, height = 6)

p_umap_annotation <- CellDimPlot(
  srt,
  group.by = "global_annotation",
  reduction = "umap",
  pt.size = 0.25,
  theme_use = "theme_blank",
  title = "GSE234129 global annotation"
)
save_plot(p_umap_annotation, "GSE234129_scop_umap_global_annotation", width = 10, height = 7)

p_umap_sample <- CellDimPlot(
  srt,
  group.by = "sample",
  reduction = "umap",
  pt.size = 0.2,
  theme_use = "theme_blank",
  title = "GSE234129 samples"
)
save_plot(p_umap_sample, "GSE234129_scop_umap_sample", width = 10, height = 7)

p_qc_violin <- FeatureStatPlot(
  srt,
  stat.by = c("n_genes_by_counts", "total_counts", "pct_counts_mt"),
  group.by = "sample",
  plot_type = "violin",
  add_box = TRUE,
  same.y.lims = FALSE,
  ncol = 1,
  title = "GSE234129 QC metrics by sample"
)
save_plot(p_qc_violin, "GSE234129_scop_qc_violin_by_sample", width = 11, height = 12)

p_qc_umap_mt <- FeatureDimPlot(
  srt,
  features = "pct_counts_mt",
  reduction = "umap",
  pt.size = 0.25,
  theme_use = "theme_blank",
  title = "GSE234129 mitochondrial percentage"
)
save_plot(p_qc_umap_mt, "GSE234129_scop_umap_pct_counts_mt", width = 8, height = 6)

p_celltype_bar <- CellStatPlot(
  srt,
  stat.by = "global_annotation",
  group.by = "sample",
  plot_type = "bar",
  stat_type = "percent",
  position = "stack",
  title = "Global annotation composition by sample",
  xlab = "Sample",
  ylab = "Percent"
)
save_plot(p_celltype_bar, "GSE234129_scop_annotation_composition_by_sample", width = 12, height = 7)

message("scop plotting complete")
