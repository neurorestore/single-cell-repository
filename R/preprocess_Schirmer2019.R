setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(Matrix)
library(Seurat)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Schirmer2019"

# read expression matrix
expr = readMM(file.path(base_dir, experiment, 'matrix.mtx'))

# add barcodes and genes
barcodes = readLines(file.path(base_dir, experiment, 'barcodes.tsv'))
genes = readLines(file.path(base_dir, experiment, 'genes.tsv'))
colnames(expr) = barcodes
rownames(expr) = gsub("^.*\t", "", genes) ## symbol only

# read metadata
meta = read.delim(file.path(base_dir, experiment, 'meta.tsv')) %>%
  # rename variables
  dplyr::rename(replicate = sample,
                label = diagnosis) %>%
  # set up global cell types
  mutate(global_cell_type = ifelse(startsWith(cell_type, "EN-"), "EN",
                                   cell_type),
         global_cell_type = ifelse(startsWith(cell_type, "IN-"), "IN",
                                   global_cell_type),
         global_cell_type = ifelse(startsWith(cell_type, "OL-"), "OL",
                                   global_cell_type)) %>%
  # barcodes as row names
  column_to_rownames('cell')

# create Seurat object
sc = CreateSeuratObject(expr, min.cells = 3, min.features = 0,
  meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
