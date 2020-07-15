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
experiment = "Kang2018"

# set up files to load
files = c("GSM2560248_2.1.mtx.gz", "GSM2560249_2.2.mtx.gz")
barcodes = c("GSM2560248_barcodes.tsv.gz", "GSM2560249_barcodes.tsv.gz")

# files are in matrix format so we can load them directly in as sparse
mats = list()
for (i in 1:length(files)) {
  mat = readMM(file.path(base_dir, experiment, files[i]))
  bars = read.delim(file.path(base_dir, experiment, barcodes[i]), header = F)
  genes = read.delim(file.path(base_dir, experiment,
                               "GSE96583_batch2.genes.tsv.gz"),
                     header = F)
  # avoid issues of same barcode across libraries
  colnames(mat) = gsub("-.*", paste0("-", i), bars$V1)
  # colnames(mat) = paste0(bars$V1, "_", i)
  rownames(mat) = genes$V2
  mats[i] = mat
}

# get merged matrix
mat = do.call(cbind, mats)

# load in meta data
meta = read.delim(file.path(base_dir, experiment,
                            "GSE96583_batch2.total.tsne.df.tsv.gz"))

# add suffix to barcodes
rows = ifelse(meta$stim == 'ctrl',
              gsub("-.*", paste0("-", 1), rownames(meta)),
              gsub("-.*", paste0("-", 2), rownames(meta)))

# flag replicate, label, and cell type
meta %<>%
  mutate(replicate = paste0('patient_', ind)) %>%
  rename(label = stim, cell_type = cell) %>%
  select(-ind) %>%
  set_rownames(rows)

# remove doublets
f1 = meta$multiplets == 'singlet'
meta = meta[f1, ]
mat = mat[, f1]

# remove unclassified cells
f2 = !is.na(meta$cell_type)
meta = meta[f2, ]
mat = mat[, f2]

# create Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
