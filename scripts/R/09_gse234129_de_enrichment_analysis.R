suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ReactomePA)
  library(msigdbr)
  library(scop)
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

marker_path <- file.path(project_dir, "results/GSE234129/tables/GSE234129_leiden05_markers.tsv")
cluster_annotation_path <- file.path(project_dir, "results/GSE234129/tables/GSE234129_cluster_annotation.tsv")
seurat_path <- file.path(project_dir, "results/GSE234129/objects/GSE234129_scop_plotting_seurat.rds")

figure_dir <- file.path(project_dir, "results/GSE234129/figures/de_enrichment")
table_dir <- file.path(project_dir, "results/GSE234129/tables/de_enrichment")
report_dir <- file.path(project_dir, "results/GSE234129/reports")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

summary_path <- file.path(table_dir, "GSE234129_de_enrichment_step_summary.tsv")
report_path <- file.path(report_dir, "GSE234129_de_enrichment_analysis_report.md")

step_results <- data.frame(step = character(), status = character(), detail = character(), stringsAsFactors = FALSE)

record_status <- function(step, status, detail = "") {
  step_results <<- rbind(step_results, data.frame(step = step, status = status, detail = detail, stringsAsFactors = FALSE))
}

run_step <- function(step, expr) {
  tryCatch(
    {
      force(expr)
      record_status(step, "completed", "")
    },
    error = function(e) {
      record_status(step, "skipped_or_failed", conditionMessage(e))
      message_time("Skipped/failed ", step, ": ", conditionMessage(e))
    }
  )
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
}

save_heatmap <- function(ht, basename, width = 10, height = 8, dpi = 600) {
  plot_obj <- if (is.list(ht) && "plot" %in% names(ht)) ht$plot else ht
  pdf_path <- file.path(figure_dir, paste0(basename, ".pdf"))
  tif_path <- file.path(figure_dir, paste0(basename, ".tif"))
  png_path <- file.path(figure_dir, paste0(basename, ".png"))
  if (inherits(plot_obj, "ggplot")) {
    ggsave(pdf_path, plot = plot_obj, width = width, height = height, device = cairo_pdf)
    ggsave(tif_path, plot = plot_obj, width = width, height = height, dpi = dpi, device = "tiff", compression = "lzw")
    ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = 180)
  } else {
    grDevices::cairo_pdf(pdf_path, width = width, height = height)
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) ComplexHeatmap::draw(plot_obj, merge_legends = TRUE) else print(plot_obj)
    grDevices::dev.off()
    grDevices::tiff(tif_path, width = width, height = height, units = "in", res = dpi, compression = "lzw")
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) ComplexHeatmap::draw(plot_obj, merge_legends = TRUE) else print(plot_obj)
    grDevices::dev.off()
    grDevices::png(png_path, width = width, height = height, units = "in", res = 180)
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) ComplexHeatmap::draw(plot_obj, merge_legends = TRUE) else print(plot_obj)
    grDevices::dev.off()
  }
  message_time("Wrote ", pdf_path)
  message_time("Wrote ", tif_path)
  message_time("Wrote ", png_path)
}

technical_gene <- function(x) {
  grepl("^MT-|^RPS|^RPL|^MTRNR|^HB[ABDEGQMZ]", x)
}

safe_neglog10 <- function(p) {
  p <- ifelse(is.na(p), NA_real_, p)
  p <- pmax(p, 1e-300)
  -log10(p)
}

normalize_data <- function(srt) {
  DefaultAssay(srt) <- "RNA"
  srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  srt
}

map_symbols <- function(genes) {
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(genes),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
}

prepare_gene_list <- function(df, split_col, min_genes = 10, max_genes = 300) {
  out <- split(df, df[[split_col]])
  out <- lapply(out, function(x) {
    x <- x[order(x$pval_adj, -x$logfoldchange), , drop = FALSE]
    unique(na.omit(x$ENTREZID))[seq_len(min(length(unique(na.omit(x$ENTREZID))), max_genes))]
  })
  out[lengths(out) >= min_genes]
}

save_enrich_result <- function(x, basename) {
  if (is.null(x)) {
    return(invisible(FALSE))
  }
  df <- as.data.frame(x)
  write_tsv(df, file.path(table_dir, paste0(basename, ".tsv")))
  nrow(df) > 0
}

plot_compare <- function(x, basename, title, width = 12, height = 8) {
  df <- as.data.frame(x)
  if (nrow(df) == 0) {
    record_status(paste0("plot_", basename), "skipped", "No enriched terms")
    return(invisible(FALSE))
  }
  p <- dotplot(x, showCategory = 8, font.size = 8, title = title) +
    theme(axis.text.y = element_text(size = 7))
  save_plot(p, basename, width = width, height = height)
  TRUE
}

message_time("Reading marker and annotation tables")
markers <- read_tsv(marker_path)
cluster_annotation <- read_tsv(cluster_annotation_path)
markers$cluster <- as.character(markers$cluster)
cluster_annotation$cluster <- as.character(cluster_annotation$cluster)

markers <- merge(
  markers,
  cluster_annotation[, c("cluster", "n_cells", "dominant_existing_celltype", "classic_marker_based_label", "proposed_global_annotation")],
  by = "cluster",
  all.x = TRUE
)
markers$broad_annotation <- sub(" \\| ref:.*$", "", markers$proposed_global_annotation)
markers$is_technical <- technical_gene(markers$gene)
markers$neg_log10_padj <- safe_neglog10(markers$pval_adj)
markers$de_status <- ifelse(
  markers$pval_adj < 0.05 & markers$logfoldchange >= 0.25,
  "up",
  ifelse(markers$pval_adj < 0.05 & markers$logfoldchange <= -0.25, "down", "not_significant")
)

mapping <- map_symbols(markers$gene)
mapping <- mapping[!is.na(mapping$ENTREZID), , drop = FALSE]
mapping <- mapping[!duplicated(mapping$SYMBOL), , drop = FALSE]
markers <- merge(markers, mapping, by.x = "gene", by.y = "SYMBOL", all.x = TRUE)

background_entrez <- unique(na.omit(markers$ENTREZID))
sig_up <- markers[
  markers$de_status == "up" & !markers$is_technical & !is.na(markers$ENTREZID),
  ,
  drop = FALSE
]

top_cluster_markers <- do.call(rbind, lapply(split(sig_up, sig_up$cluster), function(x) {
  x <- x[order(x$pval_adj, -x$logfoldchange), , drop = FALSE]
  head(x, 30)
}))
top_global_markers <- do.call(rbind, lapply(split(sig_up, sig_up$proposed_global_annotation), function(x) {
  x <- x[order(x$pval_adj, -x$logfoldchange), , drop = FALSE]
  head(x, 30)
}))
top_broad_markers <- do.call(rbind, lapply(split(sig_up, sig_up$broad_annotation), function(x) {
  x <- x[order(x$pval_adj, -x$logfoldchange), , drop = FALSE]
  head(x, 50)
}))

de_summary <- aggregate(
  cbind(n_up = markers$de_status == "up", n_down = markers$de_status == "down", n_tested = !is.na(markers$pval_adj)),
  by = list(cluster = markers$cluster, proposed_global_annotation = markers$proposed_global_annotation, broad_annotation = markers$broad_annotation),
  FUN = sum
)
de_summary$n_sig_total <- de_summary$n_up + de_summary$n_down

write_tsv(markers, file.path(table_dir, "GSE234129_de_all_leiden05_markers_annotated.tsv"))
write_tsv(de_summary, file.path(table_dir, "GSE234129_de_leiden05_summary.tsv"))
write_tsv(top_cluster_markers, file.path(table_dir, "GSE234129_de_top30_up_markers_by_cluster.tsv"))
write_tsv(top_global_markers, file.path(table_dir, "GSE234129_de_top30_up_markers_by_global_annotation.tsv"))
write_tsv(top_broad_markers, file.path(table_dir, "GSE234129_de_top50_up_markers_by_broad_annotation.tsv"))
write_tsv(mapping, file.path(table_dir, "GSE234129_de_symbol_to_entrez_mapping.tsv"))

run_step("volcano_faceted_clusters", {
  label_df <- do.call(rbind, lapply(split(sig_up, sig_up$cluster), function(x) {
    x <- x[order(x$pval_adj, -x$logfoldchange), , drop = FALSE]
    head(x, 3)
  }))
  p <- ggplot(markers, aes(x = logfoldchange, y = neg_log10_padj, color = de_status)) +
    geom_point(size = 0.22, alpha = 0.35) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = 2, linewidth = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.2) +
    geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 2,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.15
    ) +
    facet_wrap(~ cluster, scales = "free_y", ncol = 5) +
    scale_color_manual(values = c(up = "#C0392B", down = "#2878B5", not_significant = "grey75")) +
    labs(
      title = "GSE234129 Leiden 0.5 differential marker volcano panels",
      x = "log fold change",
      y = "-log10 adjusted p-value",
      color = "DE status"
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "bottom")
  save_plot(p, "GSE234129_de_volcano_faceted_leiden05", width = 14, height = 10)
})

run_step("de_summary_barplot", {
  p <- ggplot(de_summary, aes(x = reorder(cluster, as.numeric(cluster)), y = n_sig_total, fill = broad_annotation)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = n_sig_total), vjust = -0.25, size = 2.5) +
    labs(
      title = "Significant marker counts by Leiden 0.5 cluster",
      x = "Leiden 0.5 cluster",
      y = "Significant markers (padj < 0.05, |logFC| >= 0.25)",
      fill = "Broad annotation"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom")
  save_plot(p, "GSE234129_de_significant_marker_counts_by_cluster", width = 11, height = 6)
})

run_step("de_marker_heatmap", {
  srt <- readRDS(seurat_path)
  srt <- normalize_data(srt)
  features <- unique(head(top_global_markers$gene, 80))
  features <- features[features %in% rownames(srt)]
  ht <- GroupHeatmap(
    srt,
    features = features,
    group.by = "global_annotation",
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
  save_heatmap(ht, "GSE234129_de_top_global_marker_group_heatmap", width = 13, height = 10)
})

cluster_gene_list <- prepare_gene_list(sig_up, "proposed_global_annotation", min_genes = 10, max_genes = 250)
broad_gene_list <- prepare_gene_list(sig_up, "broad_annotation", min_genes = 10, max_genes = 350)

run_step("enrich_go_bp_global_annotation", {
  ego <- compareCluster(
    geneCluster = cluster_gene_list,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    universe = background_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  if (save_enrich_result(ego, "GSE234129_enrichment_GO_BP_by_global_annotation")) {
    plot_compare(ego, "GSE234129_enrichment_GO_BP_dotplot_global_annotation", "GO Biological Process enrichment by global annotation", width = 14, height = 10)
  }
})

run_step("enrich_kegg_broad_annotation", {
  ekegg <- compareCluster(
    geneCluster = broad_gene_list,
    fun = "enrichKEGG",
    organism = "hsa",
    universe = background_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20
  )
  if (save_enrich_result(ekegg, "GSE234129_enrichment_KEGG_by_broad_annotation")) {
    plot_compare(ekegg, "GSE234129_enrichment_KEGG_dotplot_broad_annotation", "KEGG enrichment by broad annotation", width = 11, height = 7)
  }
})

run_step("enrich_reactome_broad_annotation", {
  er <- compareCluster(
    geneCluster = broad_gene_list,
    fun = "enrichPathway",
    organism = "human",
    universe = background_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  if (save_enrich_result(er, "GSE234129_enrichment_Reactome_by_broad_annotation")) {
    plot_compare(er, "GSE234129_enrichment_Reactome_dotplot_broad_annotation", "Reactome enrichment by broad annotation", width = 11, height = 7)
  }
})

run_step("enrich_hallmark_broad_annotation", {
  msig_h <- tryCatch(
    msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr(species = "Homo sapiens", category = "H")
  )
  term_col <- if ("gs_name" %in% colnames(msig_h)) "gs_name" else "gs_exact_source"
  gene_col <- intersect(c("entrez_gene", "ncbi_gene", "db_ncbi_gene"), colnames(msig_h))[1]
  if (is.na(gene_col)) {
    stop("No Entrez/NCBI gene column found in msigdbr output")
  }
  hallmark <- msig_h[, c(term_col, gene_col), drop = FALSE]
  colnames(hallmark) <- c("term", "gene")
  hallmark <- hallmark[!is.na(hallmark$gene), , drop = FALSE]
  hallmark$gene <- as.character(hallmark$gene)
  eh <- compareCluster(
    geneCluster = broad_gene_list,
    fun = "enricher",
    TERM2GENE = hallmark,
    universe = background_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.25
  )
  if (save_enrich_result(eh, "GSE234129_enrichment_Hallmark_by_broad_annotation")) {
    plot_compare(eh, "GSE234129_enrichment_Hallmark_dotplot_broad_annotation", "MSigDB Hallmark enrichment by broad annotation", width = 11, height = 7)
  }
})

run_step("enrichment_network_go_bp_broad", {
  ego_broad <- compareCluster(
    geneCluster = broad_gene_list,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    universe = background_entrez,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.20,
    readable = TRUE
  )
  save_enrich_result(ego_broad, "GSE234129_enrichment_GO_BP_by_broad_annotation")
  if (nrow(as.data.frame(ego_broad)) > 0) {
    ego_broad <- pairwise_termsim(ego_broad)
    p1 <- emapplot(ego_broad, showCategory = 25) +
      ggtitle("GO BP enrichment map by broad annotation")
    save_plot(p1, "GSE234129_enrichment_GO_BP_emap_broad_annotation", width = 11, height = 9)
  }
})

write_tsv(step_results, summary_path)

figure_files <- list.files(figure_dir, pattern = "\\.(pdf|tif)$", full.names = TRUE)
table_files <- list.files(table_dir, pattern = "\\.tsv$", full.names = TRUE)

report_lines <- c(
  "# GSE234129 Differential Expression and Enrichment Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  paste0("- Marker table: `", marker_path, "`"),
  paste0("- Cluster annotation: `", cluster_annotation_path, "`"),
  paste0("- Seurat plotting object: `", seurat_path, "`"),
  "",
  "## Differential Expression",
  "",
  "- Source DE statistics are the existing Scanpy Wilcoxon Leiden 0.5 marker results.",
  "- Significant markers are defined as `pval_adj < 0.05` and `abs(logfoldchange) >= 0.25`.",
  "- Technical genes matching mitochondrial, ribosomal, or hemoglobin prefixes are excluded from enrichment gene sets.",
  paste0("- Tested marker rows: ", nrow(markers)),
  paste0("- Unique tested genes: ", length(unique(markers$gene))),
  paste0("- Mapped Entrez background genes: ", length(background_entrez)),
  paste0("- Significant up marker rows used for ORA: ", nrow(sig_up)),
  "",
  "## Enrichment",
  "",
  "- GO Biological Process enrichment was run per `global_annotation` and per broad annotation.",
  "- KEGG, Reactome, and MSigDB Hallmark ORA were run per broad annotation.",
  "- ORA background is all tested genes from the marker table with Entrez mappings.",
  "",
  "## Outputs",
  "",
  paste0("- Figure directory: `", figure_dir, "`"),
  paste0("- Table directory: `", table_dir, "`"),
  paste0("- Step summary: `", summary_path, "`"),
  paste0("- Figure files written (PDF/TIF): ", length(figure_files)),
  paste0("- Table files written: ", length(table_files)),
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
message_time("DE and enrichment analysis complete")
