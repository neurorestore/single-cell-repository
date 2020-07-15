setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(Matrix)
library(Seurat)
library(data.table)

# define base directory for raw files to preprocess
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Aztekin2019"

# read expression matrix
mat = readMM(file.path(base_dir, experiment, 'ArrayExpress/countsMatrix.mtx'))

# add barcodes and genes
barcodes = read.csv(file.path(base_dir, experiment, 'ArrayExpress/cells.csv'),
                    header = F)[[1]]
genes = read.csv(file.path(base_dir, experiment, 'ArrayExpress/genes.csv'),
                 header = F)[[1]]
colnames(mat) = barcodes
rownames(mat) = genes

# read metadata
meta = read.csv(file.path(base_dir, experiment, 'ArrayExpress/meta.csv'))
# read sample metadata
samples = read.csv(file.path(base_dir, experiment, 'ArrayExpress/labels.csv')) %>%
  set_colnames(tolower(colnames(.)))
meta %<>% left_join(samples, by = 'sample')

# tag replicate and label
meta %<>%
  dplyr::rename(replicate = sample,
                label = condition,
                cell_type = cluster) %>%
  # barcodes as row names
  column_to_rownames('cell')

# create Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
  meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
