setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(Matrix)
library(Seurat)
library(readxl)
library(loomR)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Davie2018"

# read in loom data
loom = connect(file.path(base_dir, experiment,
                         "Aerts_Fly_AdultBrain_Filtered_57k.loom"),
                       skip.validate = T)
mat = Matrix(loom[['matrix']][,], sparse = T)
genes = loom[["row_attrs/Gene"]][]
cells = loom[["col_attrs/CellID"]][]
dimnames(mat) = list(cells, genes)
mat %<>% t()

# load in cluster/cell_type annotations
anno = read.delim(file.path(base_dir, experiment, "Davie_TableS2.txt")) %>%
  set_colnames(c("cluster", "cell_type"))
# construct meta data
meta = data.frame(
  barcode = cells,
  cluster = loom[["col_attrs/Clusterings"]][]$`0`,
  label = loom[["col_attrs/Age"]][],
  replicate = loom[["col_attrs/Replicate"]][],
  gender = loom[["col_attrs/Gender"]][],
  genotype = loom[["col_attrs/Genotype"]][]
) %>%
  left_join(anno) %>%
  # adjust unannotated cell types to their original cluster
  mutate(cell_type = ifelse(cell_type == 'Unannotated', cluster, cell_type)) %>%
  column_to_rownames(var = 'barcode')

# subset matrix
mat = mat[, rownames(meta)]

# now back into Seurat
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)
# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
