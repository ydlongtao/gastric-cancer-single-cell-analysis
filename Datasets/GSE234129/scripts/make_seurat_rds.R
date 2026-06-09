suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(FALSE)
file_arg <- args[grepl("^--file=", args)]
script_path <- sub("^--file=", "", file_arg[[1]])
dataset_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
raw_dir <- file.path(dataset_dir, "raw")
out_dir <- file.path(dataset_dir, "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(raw_dir, "GSE234129_count_matrix.mtx.gz")
barcodes_path <- file.path(raw_dir, "GSE234129_barcodes.tsv.gz")
features_path <- file.path(raw_dir, "GSE234129_features.tsv.gz")
meta_path <- file.path(raw_dir, "GSE234129_meta.tsv.gz")
out_path <- file.path(out_dir, "GSE234129_seurat.rds")

counts <- readMM(matrix_path)
features <- read.delim(features_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
barcodes <- read.delim(barcodes_path, header = FALSE, stringsAsFactors = FALSE)[[1]]
metadata <- read.delim(meta_path, header = TRUE, stringsAsFactors = FALSE)

if (nrow(counts) != length(features)) {
  stop(sprintf("Feature count mismatch: matrix=%d, features=%d", nrow(counts), length(features)))
}
if (ncol(counts) != length(barcodes)) {
  stop(sprintf("Barcode count mismatch: matrix=%d, barcodes=%d", ncol(counts), length(barcodes)))
}
if (!"cell_barcodes" %in% colnames(metadata)) {
  stop("Metadata file does not contain a 'cell_barcodes' column")
}

rownames(counts) <- make.unique(features)
colnames(counts) <- barcodes
rownames(metadata) <- metadata$cell_barcodes

missing_cells <- setdiff(colnames(counts), rownames(metadata))
if (length(missing_cells) > 0) {
  stop(sprintf("Metadata missing %d cells; first missing cell: %s", length(missing_cells), missing_cells[[1]]))
}

metadata <- metadata[colnames(counts), , drop = FALSE]
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "GSE234129",
  meta.data = metadata
)
seurat_obj[["RNA"]]@misc$source <- "GEO GSE234129"

saveRDS(seurat_obj, out_path)
cat(out_path, "\n")
print(seurat_obj)
