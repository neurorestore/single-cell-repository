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
experiment = "Giladi2020"

# load in expression matrices
mat = list.files(file.path(base_dir, experiment), pattern = "AB",
  full.names = T) %>%
  map(fread) %>%
  map(column_to_rownames, var = 'V1') %>%
  map(as.matrix) %>%
  map(Matrix, sparse = T)
mat = do.call(cbind, mat)

# load in meta data files
anno = fread(file.path(base_dir, experiment, "GSE144317_sc_annotations.txt.gz"))
meta = fread(file.path(base_dir, experiment, "GSE144317_metadata.txt.gz"),
  skip = 15) %>%
  dplyr::rename(cell = well) %>%
  left_join(anno) %>%
  # remove empty wells
  filter(Number_of_cells == 1) %>%
  # remove cells with a treatment annotation
  filter(treatment != '') %>%
  # remove cells without an annotation
  filter(!is.na(annotation)) %>%
  # adjust labels
  dplyr::rename(cell_type = annotation, label = EAE.state) %>%
  # clean up some unneeded columns
  dplyr::select(-tissue, -treatment, -Number_of_cells, -filtered_in) %>%
  # set rownames
  column_to_rownames(var = 'cell')

# select relevant matrix
mat %<>% extract(, rownames(meta))

# create Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                             meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
