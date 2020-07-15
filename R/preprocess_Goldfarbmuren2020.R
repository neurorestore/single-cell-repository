setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(Matrix)
library(Seurat)
library(magrittr)
library(data.table)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Goldfarbmuren2020"

# list matrix files
mat_files = list.files(file.path(base_dir, experiment),
  pattern = "*_expression_matrix.txt$", full.names = T)

meta = fread(file.path(base_dir, experiment,
  "GSE134174_Processed_invivo_metadata.txt.gz"))

mats = mat_files %>%
  map(~ {
    id = gsub("_.*", "", basename(.))
    message(id)
    barcodes = meta %>% filter(Donor == id) %>% pull(Cell)
    barcodes = gsub("_.*", "", barcodes)
    tmp = fread(.) %>%
      column_to_rownames(var = 'V1') %>%
      dplyr::select(barcodes[barcodes %in% colnames(.)]) %>%
      as.matrix() %>%
      Matrix(sparse = T)
    colnames(tmp) = paste0(colnames(tmp), "_", id)
    return(tmp)
    })

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

# format meta data
meta %<>%
  dplyr::rename(replicate = Donor, label = Smoke_status,
    global_cell_type = cluster_ident, cell_type = subcluster_ident) %>%
  mutate(Cell = gsub("_.*", "", Cell)) %>%
  mutate(Cell = paste0(Cell, "_", replicate)) %>%
  filter(Cell %in% colnames(mat))

# remove excluded patients
meta %<>%
  filter(!grepl("excluded", label)) %>%
  droplevels() %>%
  column_to_rownames(var = 'Cell')
mat %<>% extract(, rownames(meta))

# create a new Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0, meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
