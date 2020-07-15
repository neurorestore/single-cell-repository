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
experiment = "Madissoon2020"

# first load in oesophagus_ts
eso = readRDS(file.path(base_dir, experiment, "oesophagus_ts.rds"))
# lung
lung = readRDS(file.path(base_dir, experiment, "lung_ts.rds"))
# spleen
spleen = readRDS(file.path(base_dir, experiment, "spleen_ts.rds"))

# for each object, re-format the meta data
eso@meta.data %<>%
  dplyr::rename(
    replicate = Donor,
    label = Time,
    cell_type = Celltypes
  ) %>%
  mutate(cell_type = paste0("eso_", cell_type)) %>%
  set_rownames(colnames(eso))

lung@meta.data %<>%
  dplyr::rename(
    replicate = Donor,
    label = Time,
    cell_type = Celltypes
  ) %>%
  mutate(cell_type = paste0("lung_", cell_type)) %>%
  set_rownames(colnames(lung))

spleen@meta.data %<>%
  dplyr::rename(
    replicate = Donor,
    label = Time,
    cell_type = Celltypes
  ) %>%
  mutate(cell_type = paste0("spleen_", cell_type)) %>%
  set_rownames(colnames(spleen))

# reduce the size of the objects by regenerating them
sc_eso = CreateSeuratObject(
    GetAssayData(eso, slot = 'counts'),
    min.cells = 3, min.features = 0,
    meta.data = eso@meta.data
  )
sc_lung = CreateSeuratObject(
    GetAssayData(lung, slot = 'counts'),
    min.cells = 3, min.features = 0,
    meta.data = lung@meta.data
  )
sc_spleen = CreateSeuratObject(
    GetAssayData(spleen, slot = 'counts'),
    min.cells = 3, min.features = 0,
    meta.data = spleen@meta.data
  )

# merge
sc = merge(sc_eso, c(sc_lung, sc_spleen))

# convert matrix to integers
sc@assays$RNA@counts@x = round(sc@assays$RNA@counts@x)

# save object
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
