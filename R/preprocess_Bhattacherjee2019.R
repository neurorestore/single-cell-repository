setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(Matrix)
library(Seurat)
library(data.table)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Bhattacherjee2019"

# read expression matrix
mat = fread(file.path(base_dir, experiment,
                      "GSE124952_expression_matrix.csv.gz")) %>%
  column_to_rownames(var = "V1") %>%
  as.matrix()

# read metadata
meta = read.csv(file.path(base_dir, experiment,
                          "GSE124952_meta_data.csv.gz")) %>%
  rename(age = DevStage, cell_type = L2_clusters, global_cell_type = CellType,
         label = stage, replicate = Sample) %>%
  column_to_rownames("X")

# save Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
