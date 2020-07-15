setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(Matrix)
library(Seurat)
library(magrittr)
library(readxl)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Mathys2019"

# read UMIs
dat = readMM(file.path(base_dir, experiment, "filtered_count_matrix.mtx"))
# read dimension names
genes = readLines(file.path(base_dir, experiment,
                             "filtered_gene_row_names.txt"))
columns = read.delim(file.path(base_dir, experiment,
                               "filtered_column_metadata.txt"))
rownames(dat) = genes
colnames(dat) = columns$TAG

# load remaining meta data
clinical = read.csv(file.path(base_dir, experiment,
                              "ROSMAP_Clinical_2019-05_v2.csv"))
paper = read_excel(file.path(base_dir, experiment,
                             "41586_2019_1195_MOESM5_ESM.xlsx"))
map = read.csv(file.path(base_dir, experiment,"id_mapping.csv")) %>%
  left_join(columns) %>%
  left_join(clinical) %>%
  left_join(paper %>% dplyr::select(-msex))

# setup meta data
## outcomes from paper (Fig. 2a): tangles, nft, gpath, gpath_3neocort,
## plaq_n, amyloid, cogn_global_lv
## also keep: sex, PMI
meta = map %>%
  dplyr::select(
    # cell, replicate, cell type, label
    TAG, Subject, broad.cell.type, Subcluster, pathology.group,
    # outcomes from Fig. 2a
    tangles, nft, gpath, gpath_3neocort, plaq_n, amyloid, cogn_global_lv,
    # sex and PMI
    msex, pmi,
    # other
    apoe_genotype, cogdx) %>%
  distinct() %>%
  dplyr::rename(label = pathology.group,
         replicate = Subject,
         cell_type = Subcluster) %>%
  mutate(label = ifelse(grepl("no-", label), "Control", "Alzheimers"),
         msex = ifelse(msex == 1, 'male', 'female'),
         cogdx = factor(cogdx),
         # add broad cell type with combined neurons
         combneuron.broad.cell.type = ifelse(broad.cell.type %in% c("Ex", "In"),
          'Neuron', broad.cell.type)) %>%
  column_to_rownames(var = 'TAG')

# create Seurat object
sc = CreateSeuratObject(dat, min.cells = 3, min.features = 0, meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
