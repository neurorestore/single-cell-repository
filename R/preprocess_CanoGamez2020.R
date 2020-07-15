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
experiment = "CanoGamez2020"

# load genes, barcodes, matrix
genes = fread(file.path(base_dir, experiment,
                        "NCOMMS-19-7936188_scRNAseq_genes.tsv"),
              header = F) %>% pull(V1)
barcodes = fread(file.path(base_dir, experiment,
                           "NCOMMS-19-7936188_scRNAseq_barcodes.tsv"),
                 header = F) %>% pull(V1)
mat = readMM(file.path(base_dir, experiment,
                       "NCOMMS-19-7936188_scRNAseq_raw_UMIs.mtx"))
rownames(mat) = genes
colnames(mat) = barcodes

# load in meta data
meta = fread(file.path(base_dir, experiment, 'NCOMMS-19-7936188_metadata.txt'))

# adjust meta data
meta %<>%
  dplyr::rename(cell_type = cell.type,
                replicate = donor.id,
                label = cytokine.condition) %>%
  set_rownames(colnames(mat))

# save
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
