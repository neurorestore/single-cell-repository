setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
options(datatable.fread.datatable=FALSE)
library(tidyverse)
library(data.table)
library(Matrix)
library(magrittr)
library(Seurat)

# define base directory for raw files to preprocess
base_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define GEO
experiment = "Gierahn2017"

# load matrix
mat = fread(file.path(base_dir, experiment, "GSE92495_MTB.txt.gz")) %>%
  column_to_rownames("V1")

# load metadata
meta = data.frame(
  sample = colnames(mat)) %>%
  separate(sample, c("cell_type", "label", "cell"), "_", remove = F) %>%
  mutate(cell_kind = 'Monocyte-derived macrophages') %>%
  # filter out low quality clusters
  filter(cell_type %in% c("Cluster1", "Cluster2", "Cluster3")) %>%
  column_to_rownames('sample')

# extract correct barcodes
mat %<>% extract(, rownames(meta))

# create Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
  meta.data = meta)

saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
