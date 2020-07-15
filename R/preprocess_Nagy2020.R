setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(Matrix)
library(Seurat)
library(magrittr)
library(data.table)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Nagy2020"

# load genes, barcodes, matrix
genes = read.csv(file.path(base_dir, experiment,
                        "GSE144136_GeneNames.csv.gz"))[,2]
barcodes = read.csv(file.path(base_dir, experiment,
                           "GSE144136_CellNames.csv.gz"))[,2]
mat = readMM(file.path(base_dir, experiment,
                       "GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz"))
rownames(mat) = genes
colnames(mat) = barcodes

# construct meta data
meta = data.frame(
  barcode = barcodes
) %>%
  separate(barcode, c("cell_type", "sample"), "\\.", remove = F) %>%
  separate(sample, c("replicate", "label", "batch", "barcode_strip")) %>%
  set_rownames(.$barcode)

# save
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
