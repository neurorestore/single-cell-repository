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
experiment = "Wang2020"

# load UMI count matrix
mat = fread(file.path(base_dir, experiment,
                      "GSE130664_merge_UMI_count.txt.gz")) %>%
  column_to_rownames(var = 'Gene')

# load meta data
meta = read.csv(file.path(base_dir, experiment,
                          "GSE130664_meta.csv")) %>%
  column_to_rownames(var = 'cell')
# extract annotated cells
mat = mat[, rownames(meta)]
# remove spike-ins
mat %<>% extract(!grepl("^ERCC-", rownames(.)), )

# fix meta data columns
meta %<>%
  dplyr::rename(
    replicate = individual,
    label = aging,
    cell_type = cluster
  ) %>%
  mutate(label = ifelse(label == 'O', "Old", "Young")) %>%
  set_rownames(colnames(mat))

# save
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
