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
experiment = "Wilk2020"

# load in processed Seurat object
sc = readRDS(file.path(base_dir, experiment, "blish_covid.seu.rds"))

# extract the expression matrix and metadata only
mat = GetAssayData(sc, slot = 'counts')
meta = sc@meta.data %>%
  dplyr::rename(cell_type = cell.type.fine, 
                label = Status, 
                replicate = Donor) %>%
  set_rownames(colnames(mat))

# create Seurat object
sc_new = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                             meta.data = meta)
saveRDS(sc_new, file = file.path(output_dir, paste0(experiment, ".rds")))
