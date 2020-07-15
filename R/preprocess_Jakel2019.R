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
experiment = "Jakel2019"

# read expression matrix
expr = fread(file.path(
  base_dir, experiment, 'GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz')) %>%
  column_to_rownames('V1') %>%
  as.matrix()

# read metadata
meta = read.delim(file.path(
  base_dir, experiment,'GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz')) %>%
  # rename variables
  dplyr::rename(replicate = Sample,
                label = Condition,
                cell_type = Celltypes) %>%
  # set up global cell types
  mutate(global_cell_type = gsub("[1-6]$", "", cell_type)) %>%
  # barcodes as row names
  column_to_rownames('Detected')

# create Seurat object
sc = CreateSeuratObject(expr, min.cells = 3, min.features = 0,
  meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
