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

# read raw files
accession = "Smillie2019"
mat_files = list.files(file.path(base_dir, accession),
                       pattern = ".mtx", full.names = T)
barcode_files = list.files(file.path(base_dir, accession),
                           pattern = "barcodes2.tsv", full.names = T)
genes_files = list.files(file.path(base_dir, accession),
                         pattern = "genes.tsv", full.names = T)
names = c("Epithelial", "Stromal", "Immune")
meta = fread(file.path(base_dir, accession, "all.meta2.txt"))

# merge multiple different matrices
merged_meta = list()
all_genes = list()
mats = list()
for (idx in seq_along(mat_files)) {
  mat = readMM(mat_files[idx])

  # set dimensions
  barcodes = readLines(barcode_files[idx])
  genes = readLines(genes_files[idx])
  rownames(mat) = genes
  colnames(mat) = barcodes

  # extract metadata
  meta0 = meta %>%
    filter(NAME %in% barcodes) %>%
    dplyr::rename(replicate = "Subject", label = "Health", cell_type = "Cluster",
           tissue = "Location") %>%
    mutate(tissue = names[idx]) %>%
    column_to_rownames("NAME")

  merged_meta[[idx]] = meta0
  all_genes[[idx]] = genes
  mats[[idx]] = mat
}
merged_meta %<>% do.call(rbind, .)
genes = Reduce(intersect, all_genes)
mats = map(mats, ~ extract(., genes, ))
expr = do.call(cbind, mats)

# create and save the merged object
sc = CreateSeuratObject(expr, min.cells = 3, min.features = 0,
                        meta.data = merged_meta)
saveRDS(sc, file = file.path(output_dir, paste0(accession, ".rds")))
