setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
options(datatable.fread.datatable=FALSE)
library(data.table)
library(tidyverse)
library(Matrix)
library(Seurat)
library(magrittr)

# define base directory for raw files to preprocess
base_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Mizoguchi2018"

mat = fread(file.path(base_dir, experiment,
                      "GSE109449_singlecell_rnaseq_gene_counts.tsv.gz")) %>%
  column_to_rownames(var = "ID_REF")

meta = fread(file.path(base_dir, experiment,
                      "GSE109449_singlecell_rnaseq_metadata.tsv.gz")) %>%
  mutate(cell_type = ifelse(CD34_protein == "+", "CD34+",
                                paste0("CD34-", "THY1", THY1_protein))) %>%
  mutate(cell_kind = "fibroblast") %>%
  select(sample_name, disease, age, joint, cell_kind, cell_type, donor) %>%
  dplyr::rename(replicate = "donor", label = disease) %>%
  column_to_rownames(var = "sample_name")

# create Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
  meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
