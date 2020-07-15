setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(Matrix)
library(Seurat)
library(magrittr)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "GuillaumetAdkins2017"

# load in expression matrices
mat_files = list.files(file.path(base_dir, experiment), pattern = "*.tsv.gz",
  full.names = T)
names = gsub("GSE85534_|.tsv.gz", "", basename(mat_files))

# generate complete matrix with meta data
mats = mat_files %>%
  map(read.delim) %>%
  map(column_to_rownames, var = 'X') %>%
  map(as.matrix) %>%
  map(Matrix, sparse = T) %>%
  map2(names, ~set_colnames(.x, paste0(.y, "|", colnames(.x))))

# merge into a single matrix
merge_sparse_matrices = function(matrices) {
  # adapted from https://stackoverflow.com/questions/43117608/r-binding-sparse-matrices-of-different-sizes-on-rows
  genes = map(matrices, rownames) %>%
    Reduce(union, .)
  for (matrix in matrices) {
    new_row_locs = match(rownames(matrix), genes)
    idxs = Matrix::which(matrix > 0, arr.ind = T)
    new_rows = new_row_locs[idxs[, 1]]
    cols = idxs[, 2]
    new_matrix = sparseMatrix(i = new_rows,
                              j = cols,
                              x = matrix@x,
                              dims = c(length(genes), max(cols)))
    colnames(new_matrix) = colnames(matrix)
    if (!exists("final_matrix")) {
      final_matrix = new_matrix
    }
    else {
      final_matrix = cbind2(final_matrix, new_matrix)
    }
  }
  rownames(final_matrix) = genes
  final_matrix
}
mat = merge_sparse_matrices(mats)

# construct meta data from sample name/barcodes
meta = data.frame(
  split = colnames(mat)
) %>%
  separate(split, c("meta", "barcode"), "\\|") %>%
  # adjust replicate information
  mutate(meta = gsub("_Exp", "", meta)) %>%
  separate(meta, c("cell_type", "label", "replicate"), "_") %>%
  # adjust missing replicates
  mutate(replicate = ifelse(is.na(replicate), 1, replicate)) %>%
  # set rownames
  set_rownames(.$barcode)

# adjust mat colnames to remove the meta data information
colnames(mat) = gsub(".*\\|", "", colnames(mat))

# create Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                             meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
