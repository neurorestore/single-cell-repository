setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(magrittr)
library(Matrix)
library(Seurat)
library(data.table)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Wagner2018"

# read chordin KO datasets
chd_files = list.files(file.path(base_dir, experiment), pattern = "*chd*",
                           full.names = T)
chd_mats = chd_files %>%
  extract(endsWith(., ".csv.gz") & !grepl("nm", .)) %>%
  map(fread) %>%
  map(as.data.frame) %>%
  map(~ column_to_rownames(., 'Row'))
chd_clusts = chd_files %>%
  extract(endsWith(., "_clustID.txt.gz")) %>%
  map(read.delim, header = F) %>%
  setNames(LETTERS[1:3]) %>%
  bind_rows(.id = 'batch')
# convert to a single matrix and one metadata
chd_mat = bind_cols(chd_mats) %>%
  as.matrix() %>%
  set_rownames(rownames(chd_mats[[1]]))
chd_meta = chd_clusts %>%
  set_colnames(c("batch", "cell_type")) %>%
  mutate(cell = colnames(chd_mat), genotype = "chordin KO")

# read tyrosinase KO datasets
tyr_files = list.files(file.path(base_dir, experiment), pattern = "*tyr*",
                       full.names = T)
tyr_mats = tyr_files %>%
  extract(endsWith(., ".csv.gz") & !grepl("nm", .)) %>%
  map(fread) %>%
  map(as.data.frame) %>%
  map(~ column_to_rownames(., 'Row'))
tyr_clusts = tyr_files %>%
  extract(endsWith(., "_clustID.txt.gz")) %>%
  map(read.delim, header = F) %>%
  setNames(LETTERS[1:3]) %>%
  bind_rows(.id = 'batch')
# convert to a single matrix and one metadata
tyr_mat = bind_cols(tyr_mats) %>%
  as.matrix() %>%
  set_rownames(rownames(tyr_mats[[1]]))
tyr_meta = tyr_clusts %>%
  set_colnames(c("batch", "cell_type")) %>%
  mutate(cell = colnames(tyr_mat), genotype = "tyrosinase KO")

# read WT 14hpf dataset
wt_mat = fread(file.path(base_dir, experiment, "GSM3067193_14hpf.csv.gz")) %>%
  column_to_rownames('Row') %>%
  as.matrix()
wt_meta = fread(file.path(base_dir, experiment,
                          "GSM3067193_14hpf_clustID.txt.gz"),
                header = F) %>%
  dplyr::rename(cell_type = V1) %>%
  mutate(cell = colnames(wt_mat), batch = 'D', genotype = 'WT') %>%
  extract(, colnames(tyr_meta))

# combine
mat = cbind(wt_mat, chd_mat, tyr_mat)
meta = rbind(wt_meta, chd_meta, tyr_meta) %>%
  rename(label = genotype)

# finally, map cluster indices to cluster names and flag replicate
names = fread(file.path(base_dir, experiment, "GSE112294_ClusterNames.csv.gz"))
map = with(names, setNames(ClusterName, ClusterID))
meta %<>% mutate(cell_type = map[as.character(cell_type)],
                 replicate = ifelse(label == "WT",
                                    NA,
                                    paste0(label, "_", batch)))  %>%
  column_to_rownames('cell')

# create Seurat object
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0, meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
