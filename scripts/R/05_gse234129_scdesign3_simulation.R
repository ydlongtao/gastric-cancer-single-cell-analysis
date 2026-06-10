suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(Seurat)
  library(scDesign3)
  library(ggplot2)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

command_args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  parsed <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) {
      stop("Arguments must use --name=value or --flag syntax: ", arg)
    }
    arg <- sub("^--", "", arg)
    if (grepl("=", arg, fixed = TRUE)) {
      parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
      parsed[[parts[1]]] <- paste(parts[-1], collapse = "=")
    } else {
      parsed[[arg]] <- TRUE
    }
  }
  parsed
}

args <- parse_args(command_args)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

script_arg <- commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))]
if (length(script_arg) > 0) {
  script_candidate <- sub("^--file=", "", script_arg[[1]])
  if (identical(script_candidate, "-")) {
    project_dir <- getwd()
  } else {
    script_path <- normalizePath(script_candidate, mustWork = TRUE)
    project_dir <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)
  }
} else {
  project_dir <- getwd()
}

input_seurat <- file.path(project_dir, "Datasets/GSE234129/processed/GSE234129_seurat.rds")
input_h5ad <- file.path(project_dir, "results/GSE234129/objects/GSE234129_annotated.h5ad")
python_bin <- "/Users/huangfulongtao/micromamba/envs/biomni_e1/bin/python"
export_script <- file.path(project_dir, "scripts/python/export_gse234129_scdesign3_inputs.py")

out_base <- file.path(project_dir, "results/GSE234129/scdesign3")
object_dir <- file.path(out_base, "objects")
table_dir <- file.path(out_base, "tables")
figure_dir <- file.path(out_base, "figures")
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path(table_dir, "GSE234129_scdesign3_input_metadata.tsv")
feature_meta_path <- file.path(table_dir, "GSE234129_scdesign3_input_features.tsv")
feature_set_path <- file.path(table_dir, "GSE234129_scdesign3_feature_set_hvg_markers.tsv")
top_marker_path <- file.path(project_dir, "results/GSE234129/tables/GSE234129_leiden05_top20_markers.tsv")

target_multiplier <- as.integer(args[["multiplier"]] %||% "4")
seed <- as.integer(args[["seed"]] %||% "20260609")
n_cores <- as.integer(args[["n-cores"]] %||% "2")
max_features <- as.integer(args[["max-features"]] %||% "200")
ncell_override <- args[["ncell"]]
pilot_ncell <- as.integer(args[["pilot-ncell"]] %||% "2000")
is_pilot <- isTRUE(args[["pilot"]])
prepare_only <- isTRUE(args[["prepare-only"]])
refresh_inputs <- isTRUE(args[["refresh-inputs"]])
skip_plots <- isTRUE(args[["skip-plots"]])
umap_cells_per_source <- as.integer(args[["umap-cells-per-source"]] %||% "5000")
use_bam <- !isTRUE(args[["no-bam"]])

classic_markers <- c(
  "CD3D", "CD3E", "NKG7", "GNLY", "CD8A", "IL7R",
  "MS4A1", "CD79A", "MZB1", "JCHAIN",
  "LYZ", "LST1", "CST3", "CD68", "MRC1",
  "EPCAM", "KRT8", "KRT18", "KRT19",
  "PECAM1", "VWF", "COL1A1", "DCN"
)

message_time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

save_plot <- function(plot, basename, width = 8, height = 6, dpi = 600) {
  if (skip_plots) {
    return(invisible(NULL))
  }
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

run_input_export <- function() {
  if (!file.exists(python_bin)) {
    stop("Python executable not found: ", python_bin)
  }
  if (!file.exists(export_script)) {
    stop("Input export script not found: ", export_script)
  }
  message_time("Exporting scDesign3 inputs from annotated h5ad")
  status <- system2(python_bin, export_script)
  if (!identical(status, 0L)) {
    stop("Python input export failed with status: ", status)
  }
}

sanitize_covariates <- function(metadata, covariates) {
  usable <- covariates[covariates %in% colnames(metadata)]
  usable <- usable[vapply(metadata[usable], function(x) length(unique(x[!is.na(x)])) > 1, logical(1))]
  if ("patient" %in% usable && "sample" %in% usable) {
    sample_to_patient <- tapply(metadata$patient, metadata$sample, function(x) length(unique(x)) == 1)
    if (all(sample_to_patient)) {
      message_time("sample is nested within patient; dropping patient from mu_formula to avoid collinearity")
      usable <- setdiff(usable, "patient")
    }
  }
  usable
}

qc_summary <- function(counts, metadata, source_label) {
  total_counts <- Matrix::colSums(counts)
  n_genes <- Matrix::colSums(counts > 0)
  zero_fraction <- 1 - (Matrix::colSums(counts > 0) / nrow(counts))
  data.frame(
    source = source_label,
    metric = c("total_counts", "n_genes_by_counts", "zero_fraction"),
    min = c(min(total_counts), min(n_genes), min(zero_fraction)),
    q01 = c(quantile(total_counts, 0.01), quantile(n_genes, 0.01), quantile(zero_fraction, 0.01)),
    q05 = c(quantile(total_counts, 0.05), quantile(n_genes, 0.05), quantile(zero_fraction, 0.05)),
    median = c(median(total_counts), median(n_genes), median(zero_fraction)),
    mean = c(mean(total_counts), mean(n_genes), mean(zero_fraction)),
    q95 = c(quantile(total_counts, 0.95), quantile(n_genes, 0.95), quantile(zero_fraction, 0.95)),
    q99 = c(quantile(total_counts, 0.99), quantile(n_genes, 0.99), quantile(zero_fraction, 0.99)),
    max = c(max(total_counts), max(n_genes), max(zero_fraction)),
    row.names = NULL,
    check.names = FALSE
  )
}

composition_table <- function(metadata, source_label, columns) {
  tables <- lapply(columns, function(column) {
    counts <- as.data.frame(table(metadata[[column]]), stringsAsFactors = FALSE)
    colnames(counts) <- c("level", "n_cells")
    counts$column <- column
    counts$source <- source_label
    counts$fraction <- counts$n_cells / sum(counts$n_cells)
    counts[, c("source", "column", "level", "n_cells", "fraction")]
  })
  do.call(rbind, tables)
}

marker_summary <- function(counts, metadata, markers, group_col, source_label) {
  present <- intersect(markers, rownames(counts))
  if (length(present) == 0) {
    return(data.frame())
  }
  expr <- Matrix::t(counts[present, , drop = FALSE])
  expr_df <- as.data.frame(as.matrix(expr))
  expr_df[[group_col]] <- metadata[[group_col]]
  rows <- lapply(present, function(gene) {
    aggregate(expr_df[[gene]], by = list(group = expr_df[[group_col]]), FUN = mean)
  })
  out <- do.call(rbind, Map(function(tbl, gene) {
    data.frame(source = source_label, group = tbl$group, gene = gene, mean_expression = tbl$x)
  }, rows, present))
  rownames(out) <- NULL
  out
}

build_umap_plot <- function(original_counts, original_metadata, simulated_counts, simulated_metadata, max_cells_per_source) {
  original_cells <- colnames(original_counts)
  simulated_cells <- colnames(simulated_counts)
  if (length(original_cells) > max_cells_per_source) {
    original_cells <- sample(original_cells, max_cells_per_source)
  }
  if (length(simulated_cells) > max_cells_per_source) {
    simulated_cells <- sample(simulated_cells, max_cells_per_source)
  }

  original_subset <- original_counts[, original_cells, drop = FALSE]
  simulated_subset <- simulated_counts[, simulated_cells, drop = FALSE]
  colnames(original_subset) <- paste0("original_", colnames(original_subset))
  colnames(simulated_subset) <- paste0("simulated_", colnames(simulated_subset))
  combined_counts <- cbind(original_subset, simulated_subset)

  original_meta <- original_metadata[match(original_cells, original_metadata$cell_barcode), , drop = FALSE]
  simulated_meta <- simulated_metadata[match(simulated_cells, simulated_metadata$cell_barcode), , drop = FALSE]
  combined_meta <- rbind(
    data.frame(cell_barcode = colnames(original_subset), source = "original_qc", global_annotation = original_meta$global_annotation),
    data.frame(cell_barcode = colnames(simulated_subset), source = "scdesign3_simulated", global_annotation = simulated_meta$global_annotation)
  )
  rownames(combined_meta) <- combined_meta$cell_barcode

  combined <- CreateSeuratObject(counts = combined_counts, meta.data = combined_meta, assay = "RNA", project = "GSE234129_scDesign3_validation")
  combined <- NormalizeData(combined, verbose = FALSE)
  combined <- FindVariableFeatures(combined, nfeatures = min(2000, nrow(combined)), verbose = FALSE)
  combined <- ScaleData(combined, features = VariableFeatures(combined), verbose = FALSE)
  combined <- RunPCA(combined, features = VariableFeatures(combined), npcs = 30, verbose = FALSE)
  combined <- RunUMAP(combined, dims = 1:20, verbose = FALSE)
  emb <- as.data.frame(Embeddings(combined, "umap"))
  colnames(emb)[seq_len(2)] <- c("UMAP_1", "UMAP_2")
  emb$source <- combined$source
  emb$global_annotation <- combined$global_annotation

  ggplot(emb, aes(x = UMAP_1, y = UMAP_2, color = source)) +
    geom_point(size = 0.18, alpha = 0.55) +
    theme_bw() +
    coord_equal() +
    labs(x = "UMAP 1", y = "UMAP 2", color = "Source", title = "GSE234129 original vs scDesign3 simulated")
}

build_marker_dotplot <- function(marker_summary_table) {
  if (nrow(marker_summary_table) == 0) {
    return(NULL)
  }
  marker_summary_table$group <- as.character(marker_summary_table$group)
  top_groups <- names(sort(table(marker_summary_table$group), decreasing = TRUE))[seq_len(min(12, length(unique(marker_summary_table$group))))]
  plot_data <- marker_summary_table[marker_summary_table$group %in% top_groups, , drop = FALSE]
  ggplot(plot_data, aes(x = gene, y = group, size = mean_expression, color = source)) +
    geom_point(alpha = 0.75) +
    facet_wrap(~source) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
    labs(x = "Marker", y = "Global annotation", size = "Mean expression", color = "Source", title = "Classic marker expression")
}

if (refresh_inputs || !file.exists(metadata_path) || !file.exists(feature_meta_path)) {
  run_input_export()
}

if (!file.exists(metadata_path) || !file.exists(feature_meta_path)) {
  stop("Missing exported input metadata/features. Run with --refresh-inputs.")
}
if (!file.exists(input_seurat)) {
  stop("Missing Seurat input object: ", input_seurat)
}
if (!file.exists(input_h5ad)) {
  stop("Missing annotated h5ad input: ", input_h5ad)
}

set.seed(seed)
message_time("Reading input metadata and feature metadata")
metadata <- read.delim(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
feature_meta <- read.delim(feature_meta_path, stringsAsFactors = FALSE, check.names = FALSE)
required_metadata <- c("cell_barcode", "patient", "sample", "celltype", "global_annotation")
missing_metadata <- setdiff(required_metadata, colnames(metadata))
if (length(missing_metadata) > 0) {
  stop("Input metadata is missing columns: ", paste(missing_metadata, collapse = ", "))
}

message_time("Reading Seurat counts object")
seurat_obj <- readRDS(input_seurat)
counts_all <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
common_cells <- intersect(metadata$cell_barcode, colnames(counts_all))
if (length(common_cells) == 0) {
  stop("No overlapping cells between annotated h5ad metadata and Seurat counts")
}
metadata <- metadata[match(common_cells, metadata$cell_barcode), , drop = FALSE]
counts_all <- counts_all[, metadata$cell_barcode, drop = FALSE]

hvg_genes <- feature_meta$gene[feature_meta$highly_variable %in% c(TRUE, "True", "true", "1", 1)]
marker_genes <- classic_markers
cluster_top_marker_genes <- character()
if (file.exists(top_marker_path)) {
  top_markers <- read.delim(top_marker_path, stringsAsFactors = FALSE)
  if ("gene" %in% colnames(top_markers)) {
    cluster_top_marker_genes <- unique(top_markers$gene)
    marker_genes <- unique(c(marker_genes, cluster_top_marker_genes))
  }
}
marker_genes <- unique(marker_genes[marker_genes %in% rownames(counts_all)])
if ("dispersions_norm" %in% colnames(feature_meta)) {
  hvg_meta <- feature_meta[feature_meta$gene %in% hvg_genes, c("gene", "dispersions_norm"), drop = FALSE]
  hvg_meta$dispersions_norm <- suppressWarnings(as.numeric(hvg_meta$dispersions_norm))
  hvg_meta <- hvg_meta[order(hvg_meta$dispersions_norm, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  hvg_genes_ranked <- hvg_meta$gene
} else {
  hvg_genes_ranked <- hvg_genes
}
hvg_genes_ranked <- unique(hvg_genes_ranked[hvg_genes_ranked %in% rownames(counts_all)])

feature_genes <- unique(c(marker_genes, hvg_genes_ranked))
feature_genes <- feature_genes[feature_genes %in% rownames(counts_all)]
if (!is.na(max_features) && max_features > 0 && length(feature_genes) > max_features) {
  feature_genes <- feature_genes[seq_len(max_features)]
}
if (length(feature_genes) == 0) {
  stop("No HVG+marker genes overlap the count matrix")
}

feature_source <- data.frame(
  gene = feature_genes,
  in_hvg = feature_genes %in% hvg_genes,
  in_classic_marker = feature_genes %in% classic_markers,
  in_cluster_top_marker = feature_genes %in% cluster_top_marker_genes,
  stringsAsFactors = FALSE
)
write.table(feature_source, feature_set_path, sep = "\t", quote = FALSE, row.names = FALSE)
message_time("Wrote feature set: ", feature_set_path)
message_time("Feature set size: ", length(feature_genes), " genes")

counts <- counts_all[feature_genes, , drop = FALSE]
metadata_for_sce <- metadata
rownames(metadata_for_sce) <- metadata_for_sce$cell_barcode
for (column in c("patient", "sample", "celltype", "global_annotation")) {
  metadata_for_sce[[column]] <- factor(metadata_for_sce[[column]])
}

baseline_ncell <- ncol(counts)
target_ncell <- baseline_ncell * target_multiplier
if (!is.null(ncell_override)) {
  target_ncell <- as.integer(ncell_override)
}
if (is_pilot) {
  target_ncell <- pilot_ncell
}
if (!is_pilot && target_ncell < baseline_ncell * 4) {
  stop("Formal run must simulate at least four times the QC baseline cells: ", baseline_ncell * 4)
}

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = S4Vectors::DataFrame(metadata_for_sce)
)

covariates <- sanitize_covariates(metadata_for_sce, c("global_annotation", "patient", "sample"))
if (!"global_annotation" %in% covariates) {
  stop("global_annotation is required as the primary scDesign3 celltype variable")
}
formula_terms <- paste(covariates, collapse = " + ")
mu_formula <- args[["mu-formula"]] %||% formula_terms
sigma_formula <- args[["sigma-formula"]] %||% "1"

run_label <- if (is_pilot) "pilot" else "4x"
feature_label <- paste0(length(feature_genes), "features")
output_prefix <- paste0("GSE234129_scdesign3_simulated_", run_label, "_", feature_label, "_hvg_markers")
counts_out <- file.path(object_dir, paste0(output_prefix, "_counts.rds"))
sce_out <- file.path(object_dir, paste0(output_prefix, "_sce.rds"))
seurat_out <- file.path(object_dir, paste0(output_prefix, "_seurat.rds"))
metadata_out <- file.path(table_dir, paste0(output_prefix, "_metadata.tsv"))
qc_summary_out <- file.path(table_dir, "GSE234129_scdesign3_original_vs_simulated_qc_summary.tsv")
composition_out <- file.path(table_dir, "GSE234129_scdesign3_original_vs_simulated_composition.tsv")
marker_summary_out <- file.path(table_dir, "GSE234129_scdesign3_original_vs_simulated_marker_summary.tsv")
report_out <- file.path(out_base, "GSE234129_scdesign3_simulation_report.md")

if (prepare_only) {
  message_time("Prepared inputs only; skipping scDesign3 because --prepare-only was set")
  quit(save = "no", status = 0)
}

message_time("Running scDesign3")
message_time("Baseline cells: ", baseline_ncell)
message_time("Target simulated cells: ", target_ncell)
message_time("mu_formula: ", mu_formula)

simu <- scdesign3(
  sce = sce,
  assay_use = "counts",
  celltype = "global_annotation",
  other_covariates = setdiff(covariates, "global_annotation"),
  ncell = target_ncell,
  mu_formula = mu_formula,
  sigma_formula = sigma_formula,
  family_use = "nb",
  n_cores = n_cores,
  usebam = use_bam,
  corr_formula = "1",
  copula = "gaussian",
  if_sparse = TRUE,
  DT = TRUE,
  pseudo_obs = FALSE,
  return_model = FALSE,
  simplify = TRUE,
  nonzerovar = TRUE,
  trace = TRUE
)

sim_counts <- simu$new_count
if (is.list(sim_counts) && !is.data.frame(sim_counts)) {
  sim_counts <- sim_counts[[1]]
}
if (!inherits(sim_counts, "Matrix")) {
  sim_counts <- Matrix(sim_counts, sparse = TRUE)
}
sim_counts <- round(sim_counts)
sim_counts[sim_counts < 0] <- 0
rownames(sim_counts) <- rownames(counts)

sim_metadata <- as.data.frame(simu$new_covariate, stringsAsFactors = FALSE)
sim_metadata$cell_barcode <- paste0("scDesign3_", seq_len(ncol(sim_counts)))
colnames(sim_counts) <- sim_metadata$cell_barcode
rownames(sim_metadata) <- sim_metadata$cell_barcode

sim_sce <- SingleCellExperiment(
  assays = list(counts = sim_counts),
  colData = S4Vectors::DataFrame(sim_metadata)
)

sim_seurat <- CreateSeuratObject(
  counts = sim_counts,
  meta.data = sim_metadata,
  assay = "RNA",
  project = "GSE234129_scDesign3"
)

saveRDS(sim_counts, counts_out)
saveRDS(sim_sce, sce_out)
saveRDS(sim_seurat, seurat_out)
write.table(sim_metadata, metadata_out, sep = "\t", quote = FALSE, row.names = FALSE)
message_time("Wrote ", counts_out)
message_time("Wrote ", sce_out)
message_time("Wrote ", seurat_out)
message_time("Wrote ", metadata_out)

original_summary <- qc_summary(counts, metadata_for_sce, "original_qc")
sim_summary <- qc_summary(sim_counts, sim_metadata, "scdesign3_simulated")
write.table(rbind(original_summary, sim_summary), qc_summary_out, sep = "\t", quote = FALSE, row.names = FALSE)

composition_columns <- intersect(c("patient", "sample", "celltype", "global_annotation"), colnames(sim_metadata))
composition <- rbind(
  composition_table(metadata_for_sce, "original_qc", composition_columns),
  composition_table(sim_metadata, "scdesign3_simulated", composition_columns)
)
write.table(composition, composition_out, sep = "\t", quote = FALSE, row.names = FALSE)

marker_group_col <- if ("global_annotation" %in% colnames(sim_metadata)) "global_annotation" else composition_columns[[1]]
marker_summary_table <- rbind(
  marker_summary(counts, metadata_for_sce, classic_markers, marker_group_col, "original_qc"),
  marker_summary(sim_counts, sim_metadata, classic_markers, marker_group_col, "scdesign3_simulated")
)
write.table(marker_summary_table, marker_summary_out, sep = "\t", quote = FALSE, row.names = FALSE)

qc_long <- rbind(
  data.frame(
    source = "original_qc",
    total_counts = Matrix::colSums(counts),
    n_genes_by_counts = Matrix::colSums(counts > 0),
    zero_fraction = 1 - Matrix::colSums(counts > 0) / nrow(counts)
  ),
  data.frame(
    source = "scdesign3_simulated",
    total_counts = Matrix::colSums(sim_counts),
    n_genes_by_counts = Matrix::colSums(sim_counts > 0),
    zero_fraction = 1 - Matrix::colSums(sim_counts > 0) / nrow(sim_counts)
  )
)

p_qc <- ggplot(qc_long, aes(x = source, y = n_genes_by_counts, fill = source)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  geom_boxplot(width = 0.16, outlier.size = 0.2, linewidth = 0.2) +
  facet_wrap(~source, scales = "free_x") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "Genes detected", title = "GSE234129 scDesign3 QC distribution")
save_plot(p_qc, "GSE234129_scdesign3_original_vs_simulated_qc", width = 8, height = 5)

if ("global_annotation" %in% composition_columns) {
  comp_subset <- composition[composition$column == "global_annotation", , drop = FALSE]
  p_comp <- ggplot(comp_subset, aes(x = source, y = fraction, fill = level)) +
    geom_col(width = 0.72) +
    theme_bw() +
    labs(x = NULL, y = "Fraction", fill = "Global annotation", title = "GSE234129 scDesign3 composition")
  save_plot(p_comp, "GSE234129_scdesign3_original_vs_simulated_composition", width = 10, height = 6)
}

if (!skip_plots) {
  message_time("Building sampled original-vs-simulated UMAP validation plot")
  p_umap <- build_umap_plot(counts, metadata, sim_counts, sim_metadata, umap_cells_per_source)
  save_plot(p_umap, "GSE234129_scdesign3_original_vs_simulated_umap", width = 8, height = 6)

  p_marker <- build_marker_dotplot(marker_summary_table)
  if (!is.null(p_marker)) {
    save_plot(p_marker, "GSE234129_scdesign3_original_vs_simulated_marker_dotplot", width = 12, height = 7)
  }
}

if (ncol(sim_counts) >= baseline_ncell * 4 || is_pilot) {
  validation_status <- "passed"
} else {
  validation_status <- "failed"
}

report_lines <- c(
  "# GSE234129 scDesign3 Simulation Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  "",
  paste0("- Annotated h5ad: `", input_h5ad, "`"),
  paste0("- Seurat counts: `", input_seurat, "`"),
  paste0("- QC baseline cells: ", baseline_ncell),
  paste0("- Feature strategy: HVG + classic markers + Leiden 0.5 top markers"),
  paste0("- Feature count: ", length(feature_genes)),
  paste0("- max_features: ", max_features),
  "",
  "## scDesign3 Parameters",
  "",
  paste0("- seed: ", seed),
  paste0("- ncell: ", target_ncell),
  paste0("- celltype: `global_annotation`"),
  paste0("- other_covariates: `", paste(setdiff(covariates, "global_annotation"), collapse = "`, `"), "`"),
  paste0("- mu_formula: `", mu_formula, "`"),
  paste0("- sigma_formula: `", sigma_formula, "`"),
  "- family_use: `nb`",
  paste0("- usebam: `", use_bam, "`"),
  "- copula: `gaussian`",
  "- DT: `TRUE`",
  "",
  "## Outputs",
  "",
  paste0("- `", counts_out, "`"),
  paste0("- `", sce_out, "`"),
  paste0("- `", seurat_out, "`"),
  paste0("- `", metadata_out, "`"),
  paste0("- `", qc_summary_out, "`"),
  paste0("- `", composition_out, "`"),
  "",
  "## Validation",
  "",
  paste0("- Simulated cells: ", ncol(sim_counts)),
  paste0("- Simulated genes: ", nrow(sim_counts)),
  paste0("- Metadata rows: ", nrow(sim_metadata)),
  paste0("- Status: ", validation_status),
  "",
  "## Session Info",
  "",
  "```text",
  capture.output(sessionInfo()),
  "```"
)
writeLines(report_lines, report_out)
message_time("Wrote ", report_out)
message_time("scDesign3 simulation complete")
