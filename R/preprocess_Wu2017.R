setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(Matrix)
library(Seurat)
library(magrittr)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECTS"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Wu2017"

# read dataset
dat = read.delim(file.path(base_dir, experiment,
                           "GSE103976_MeA_AllCells_DGE.txt.gz"))

# read neuronal dataset, to obtain subtypes
neu = read.delim(file.path(base_dir, experiment,
                           "GSE103976_MeA_Neurons_DGE.txt.gz"))

# construct metadata data frame from column names
cols = colnames(dat)
not_barcodes = gsub("^.*[AGTCN]_", "", cols)
cell_types = gsub("^.*_", "", not_barcodes)
condition = gsub("\\..*$", "", not_barcodes)
replicates = gsub("^.*\\.|_.*$", "", not_barcodes) %>%
  paste0(condition, "_", .)
barcodes = paste0(gsub("_.*$", "", cols), "_", condition, ".", replicates)
meta = data.frame(cell_type = cell_types,
                  replicate = replicates,
                  condition = condition,
                  barcode = barcodes,
                  names = cols) %>%
  mutate(condition = fct_recode(condition,
                                'stress' = 'IM',
                                'odor' = 'O',
                                'seizure' = 'SZ',
                                'control' = 'CTL'))

# do the same for neurons and merge
neu_cols = colnames(neu)[-1]
neu_not_barcodes = gsub("^.*[AGTCN]_", "", neu_cols)
neu_cell_types = gsub("^.*_", "", neu_not_barcodes)
neu_replicates = gsub("_.*$", "", neu_not_barcodes)
neu_replicates = gsub("\\.", "_", neu_replicates)
neu_condition = gsub("\\..*$", "", neu_not_barcodes)
neu_barcodes = paste0(gsub("_.*$", "", neu_cols), "_", neu_condition, ".",
                      neu_replicates)
neu_meta = data.frame(cell_type = neu_cell_types,
                      replicate = neu_replicates,
                      condition = neu_condition,
                      barcode = neu_barcodes,
                      names = neu_cols) %>%
  mutate(condition = fct_recode(condition,
                                'stress' = 'IM',
                                'odor' = 'O',
                                'seizure' = 'SZ',
                                'control' = 'CTL'))
meta %<>%
  left_join(neu_meta %>%
    rename(cell_subtype = cell_type) %>%
    select(cell_subtype, barcode)
)

# replace cell type with subtype for neurons
meta %<>% mutate(cell_type = ifelse(cell_type == "N" & !is.na(cell_subtype),
  cell_subtype, cell_type)) %>%
  dplyr::rename(label = condition)

# convert counts to matrix
mat = as.matrix(dat)

# create Seurat object
meta %<>%
  remove_rownames() %>%
  column_to_rownames('names')
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
