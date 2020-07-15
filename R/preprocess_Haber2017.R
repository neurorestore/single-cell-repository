setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(Matrix)
library(Seurat)
library(magrittr)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Haber2017"

# load droplet matrix
drop = read.delim(file.path(base_dir, experiment,
                            "GSE92332_SalmHelm_UMIcounts.txt.gz"))

# parse meta data for each cell from column names
drop_meta = data.frame(sample = colnames(drop)) %>%
  separate(sample, c("replicate", "barcode", "label", "cell_type"),
           sep = "_", remove = F) %>%
  column_to_rownames(var = 'sample')

# create Seurat object
sc_drop = CreateSeuratObject(drop, min.cells = 3, min.features = 0,
                             meta.data = drop_meta)
saveRDS(sc_drop, file = file.path(output_dir, paste0(experiment, ".rds")))
