setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(Matrix)
library(Seurat)
library(sctransform)
library(openxlsx)

# define base directory for raw files to preprocess
base_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define GEO accession
experiment = "McGinni2019"

# read expression matrix
expr_dataset = 'HMEC_orig'
message(".. loading matrix for: ", expr_dataset)

# read matrix
mat = readMM(file.path(base_dir, experiment,
                        paste0('GSE129578_', expr_dataset, '_matrix.mtx.gz')))

# add barcodes and genes
barcodes = readLines(file.path(base_dir, experiment,
                               paste0(expr_dataset, '_barcodes.tsv')))
genes = readLines(file.path(base_dir, experiment,
                            paste0(expr_dataset, '_genes.tsv')))
colnames(mat) = barcodes
rownames(mat) = gsub("^.*\t", "", genes) ## symbol only

# read metadata
meta = openxlsx::read.xlsx(
  file.path(base_dir, experiment, '41592_2019_433_MOESM9_ESM.xlsx'), startRow = 2) %>%
  dplyr::rename(cell_type = Composition,
                label = Condition) %>%
  remove_rownames() %>%
  column_to_rownames('CellID')

# subset expression matrix
keep = rownames(meta) %>%
  intersect(colnames(mat))
mat %<>% extract(, keep)
meta %<>% extract(keep, )

# create Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
