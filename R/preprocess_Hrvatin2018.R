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
experiment = "Hrvatin2018"

# read count matrix
count_files = list.files(file.path(base_dir, experiment),
                         pattern = 'counts', full.names = T) %>%
  extract(grepl("0h|1h|4h", .))
dats = map(count_files, fread)

# read metadata
meta = fread(file.path(base_dir, experiment,
                       "GSE102827_cell_type_assignments.csv.gz"),
             header = T)

# match up barcodes between metadata and count matrices
dats0 = list()
for (idx in seq_along(count_files)) {
  dat = dats[[idx]]
  count_file = count_files[idx]
  sample_name = basename(count_file) %>%
    strsplit('_') %>%
    unlist() %>%
    extract(seq(2, length(.) - 1)) %>%
    paste0(collapse = '_')
  barcodes = paste0(sample_name, '_', colnames(dat))
  colnames(dat) = barcodes
  dats0[[idx]] = dat
}
meta %<>% mutate(barcode = paste0(sample, "_", gsub("^.*_", "", V1)))
expr = bind_cols(dats0)

# fix the rest of the metadata
meta %<>%
  filter(!is.na(maintype)) %>%
  mutate(cell_type = ifelse(is.na(subtype), celltype, subtype)) %>%
  rename(label = stim) %>%
  mutate(replicate = gsub("^B*", "", sample)) %>%
  mutate(replicate = gsub("^._", "", replicate)) %>%
  mutate(replicate = gsub("_.*", "", replicate)) %>%
  mutate(replicate = paste0("animal_", replicate)) %>%
  column_to_rownames("barcode")

# filter barcodes without metadata
genes = expr$B1_1_0h_gene
expr %<>%
  dplyr::select(rownames(meta)) %>%
  as.matrix() %>%
  set_rownames(genes)

# save Seurat object
sc = CreateSeuratObject(expr, min.cells = 3, min.features = 0,
                        meta.data = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
