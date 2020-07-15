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
experiment = "Angelidis2019"

# read the Seurat object with raw counts (from GEO upload)
input_dir = file.path(base_dir, experiment)
load(file.path(input_dir, "GSE124872_raw_counts_single_cell.RData"),
     verbose = T)
expr = raw_counts ## rename

# match the counts up to metadata from GitHub
load(file.path(input_dir, 'data', 'SeuratObject.RData'), verbose = T)
sc = seu.ica
meta = sc@meta.data %>%
  rownames_to_column('barcode')

# we need to construct a map between the two sets of barcodes,
# to exclude the possibility they are not ordered differently
bc1 = colnames(expr) %>% gsub("^.*:", "", .)
bc2 = meta$barcode %>% gsub("^.*:", "", .)
all(bc1 == bc2) ## TRUE

# having confirmed this, we can standardize the metadata ...
meta %<>%
  dplyr::rename(replicate = identifier,
                label = grouping,
                cell_type = celltype) %>%
  droplevels() %>%
  mutate(label = as.character(label)) %>%
  dplyr::select(replicate, label, cluster, cell_type)
rownames(meta) = colnames(expr)

# remove two cells with missing types
keep = !is.na(meta$cell_type)
meta %<>% extract(keep, )
expr %<>% extract(, keep)

# reconstruct Seurat object
sc = CreateSeuratObject(expr, min.cells = 3, min.features = 0,
  meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
