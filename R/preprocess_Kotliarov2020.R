## NOTE: Run this with r/3.5.0 and Seurat V2
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
experiment = "Kotliarov2020"

# load in seurat object provided by authors online
sc = readRDS(file.path(base_dir, experiment,
                       "H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds"))

# load in cluster annotations
anno = read.delim(file.path(base_dir, experiment,
                            'clustree_node_labels_withCellTypeLabels.txt')) %>%
  filter(clustering == 'p3_dist_3') %>%
  dplyr::rename(K3 = label)

# add this to the meta data
sc@meta.data %<>%
  left_join(anno) %>%
  set_rownames(colnames(sc@data))

# regenerate Seurat object to reduce size
expr = sc@raw.data[, colnames(sc@assay$CITE@data)]
sc_new = CreateSeuratObject(
  expr, min.cells = 3, min.features = 0,
  meta = data.frame(
    batch = sc@meta.data$batch,
    replicate = sc@meta.data$sampleid,
    label = gsub("d0 ", "", sc@meta.data$adjmfc.time),
    cell_type = sc@meta.data$Cell.Type.label
  ) %>%
    set_rownames(colnames(sc@data))
)
saveRDS(sc_new, file = file.path(output_dir, paste0(experiment, ".rds")))

## quit the R session, load v3 back in
sc = readRDS(file = file.path(output_dir, paste0(experiment, ".rds")))
sc_new = UpdateSeuratObject(sc)
# re-save it
saveRDS(sc_new, file = file.path(output_dir, paste0(experiment, ".rds")))
