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
experiment = "Huang2020"

# read expression matrix
mat = fread(file.path(base_dir, experiment,
                      "Pediatric_IBD_colitis_scRNA_countMatrix.txt")) %>%
  column_to_rownames(var = "V1")

# read metadata
meta = fread(file.path(base_dir, experiment,
                          "Pediatric_IBD_colitis_scRNA_annotated.txt")) %>%
  rename(replicate = individual, label = group, global_cell_type = major,
         cell_type = subset) %>%
  column_to_rownames("V1")

# save Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
