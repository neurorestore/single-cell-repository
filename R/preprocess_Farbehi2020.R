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
experiment = "Farbehi2020"

# load in droplet matrices for FACS and flow through cells
mat1 = fread(file.path(base_dir, experiment, "GFP_ShamVsMI_days3_7.txt")) %>%
  column_to_rownames(var = 'V1') %>%
  as.matrix() %>%
  Matrix(sparse = T)
mat2 = fread(file.path(base_dir, experiment, "TIP_ShamVsMI_days3_7.txt")) %>%
  column_to_rownames(var = 'V1') %>%
  as.matrix() %>%
  Matrix(sparse = T)

# load in meta data
meta1 = fread(file.path(base_dir, experiment, "GFP_tSNE_cluster_ID_table.txt")) %>%
  dplyr::select(cell_barcode, cluster, experiment) %>%
  mutate(type = "GFP") %>%
  dplyr::rename(cell_type = cluster, label = experiment) %>%
  mutate(label = gsub("day", "Day", label)) %>%
  column_to_rownames(var = 'cell_barcode')
meta2 = fread(file.path(base_dir, experiment, "TIP_tSNE_cluster_ID_table.txt")) %>%
  dplyr::select(cell_barcode, cluster, experiment) %>%
  mutate(type = "TIP") %>%
  dplyr::rename(cell_type = cluster, label = experiment) %>%
  column_to_rownames(var = 'cell_barcode')

mat1 %<>% extract(, rownames(meta1))
mat2 %<>% extract(, rownames(meta2))

# create Seurat objects
sc1 = CreateSeuratObject(mat1, min.cells = 3, min.features = 0, meta.data = meta1)
sc2 = CreateSeuratObject(mat2, min.cells = 3, min.features = 0, meta.data = meta2)

# merge Seurat objects and save
sc_final = merge(sc1, sc2)
saveRDS(sc_final, file = file.path(output_dir, paste0(experiment, ".rds")))
