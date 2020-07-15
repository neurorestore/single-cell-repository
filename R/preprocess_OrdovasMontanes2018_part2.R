# Part 2: pick up preprocessing where the Supplementary Data script left off
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
experiment = "OrdovasMontanes2018"

data = readRDS(file.path(base_dir, experiment, "part_1",
                         "OrdovasMontanes2018_data.rds"))
meta = readRDS(file.path(base_dir, experiment, "part_1",
                         "OrdovasMontanes2018_meta.rds")) %>%
  rownames_to_column(var = "names") %>%
  rename(label = "polyp", cell_type = subset) %>%
  mutate(label = ifelse(label == "YES", "polyp", "healthy")) %>%
  mutate(replicate = gsub("TOT", "", orig.ident)) %>%
  column_to_rownames(var = "names")

# create Seurat object
sc = CreateSeuratObject(data, min.cells = 3, min.features = 0,
                        meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
