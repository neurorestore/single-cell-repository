setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(Matrix)
library(Seurat)
library(data.table)
library(magrittr)
library(reticulate)
np = import('numpy')

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Schafflick2020"

# construct PBMC matrix
## only need one genes file (they are all the same)
gene_files = list.files(file.path(base_dir, experiment),
                        pattern = "*_genes.tsv.gz", full.names = T)
genes = fread(gene_files[1], header = F) %>% pull(V1)

# define order of patients to ensure meta data lines up
patients = c('MS19270', 'MS49131','MS71658','MS60249','MS74594',
             'PTC32190', 'PTC41540','PTC85037','PST83775','PST95809')
gsm = c('GSM4104134', 'GSM4104136','GSM4104135','GSM4104137','GSM4104138',
        'GSM4104140', 'GSM4104142','GSM4104143','GSM4104139','GSM4104141')

## get barcodes
barcode_files = file.path(base_dir, experiment,
                          paste0(gsm, "_", patients,
                                 "_PBMCs_GRCh38_barcodes.tsv.gz"))
barcodes = barcode_files %>%
  map(fread, header = F) %>%
  setNames(basename(barcode_files)) %>%
  bind_rows(.id = 'file') %>%
  mutate(barcode = paste0(V1, "_", file)) %>%
  pull(barcode)

## get matrices
mat_files = file.path(base_dir, experiment,
                      paste0(gsm, "_", patients,
                             "_PBMCs_GRCh38_matrix.mtx.gz"))
mats = map(mat_files, readMM)
mat = do.call(cbind, mats)
colnames(mat) = barcodes
rownames(mat) = genes

# save Seurat object
## PBMCs
sc_pbmc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                             meta.data = pbmc_meta)
saveRDS(sc_pbmc, file = file.path(output_dir, paste0(experiment, ".rds")))
