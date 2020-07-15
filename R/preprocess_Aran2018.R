setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(Matrix)
library(Seurat)
library(magrittr)
library(data.table)


# define base directory for raw files to preprocess
base_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Aran2018"

# load singler object obtained from https://github.com/dviraran/SingleR/blob/master/manuscript_figures/FiguresData/GSE111664.RData
load(file.path(base_dir, experiment, 'GSE111664.RData'))
anno = data.frame(cell_type = singler$singler[[1]]$SingleR.clusters$labels[,],
                  clusters = names(singler$singler[[1]]$SingleR.clusters$labels[,]))
meta = singler$seurat@meta.data %>%
  dplyr::rename(clusters = "res.0.8", replicate = "orig.ident") %>%
  mutate(label = ifelse(grepl('Bleomycin', replicate),
                        'fibrosis', 'control')) %>%
  left_join(anno)
rownames(meta) = rownames(singler$seurat@meta.data)

# set barcodes to match with GEO deposit
meta %<>% mutate(barcode = paste0(replicate, "-", rownames(.)))

# load in GEO matrices, bind, convert to sparse
names = c(paste0("Control_", c(1, seq(1,6))), paste0("Bleomycin_", seq(1,3)))
mats = list.files(file.path(base_dir, experiment), pattern = "*.txt.gz",
  full.names = T) %>%
  map(fread) %>%
  map(column_to_rownames, var = 'GENE') %>%
  map2(names, ~ set_colnames(.x, paste0(.y, "-", colnames(.x)))) %>%
  map(as.matrix) %>%
  map(Matrix, sparse = T)

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

# barcodes to keep
keep = meta$barcode[meta$barcode %in% colnames(mat)]
meta %<>% filter(barcode %in% keep)
mat %<>% extract(, keep)

sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
  meta.data = meta)

saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
