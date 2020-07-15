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
experiment = "Hu2017"

# read processed matrix
## from https://github.com/wulabupenn/Hu_MolCell_2017/blob/master/MolCell_Fig4.md
dat = readRDS(file.path(base_dir, experiment, "data", "scNuc2.rds"))

# move out of Seurat
mat = dat@raw.data
meta = dat@data.info %>%
  rownames_to_column('names')

# keep only PTZ vs. saline comparison
meta %<>% filter(type %in% c("PTZ", "saline"))
keep = match(meta$cell, colnames(mat))
mat %<>% extract(, keep)

# keep clustering at relevant resolution (cf. Github)
meta %<>%
  dplyr::select(cell, type, ID, res.comb2, names) %>%
  set_colnames(c("barcode", "label", "replicate", "cell_type", "names")) %>%
  mutate(replicate = as.character(replicate)) %>%
  column_to_rownames("names")

# now back into Seurat
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
