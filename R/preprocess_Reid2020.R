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
experiment = "Reid2020"

#### First do P.Bergei
# load the count data
mat = read_excel(file.path(base_dir, experiment, "elife-33105-supp6-v1.xlsx"),
                 sheet = 1)

# load meta data for this experiment
meta = read_excel(file.path(base_dir, experiment, "elife-33105-supp5-v1.xlsx"),
                  sheet = 4) %>%
  # remove cells that did not pass the quality filter
  filter(pass_filter) %>%
  dplyr::rename(label = consensus) %>%
  mutate(cell_type = 'p.berghei') %>%
  column_to_rownames(var = 'sample_id')

# subset mat
mat0 = mat[, c('id', rownames(meta))]
mat0 %<>% column_to_rownames(var = 'id')

sc = CreateSeuratObject(mat0, meta = meta)
saveRDS(sc, file = file.path(output_dir, paste0(experiment, "_berghei.rds")))

#### Next do P.falciparum asexual
# load the count data
mat2 = read_excel(file.path(base_dir, experiment, "elife-33105-supp6-v1.xlsx"),
                  sheet = 2)

# load meta data for this experiment
meta2 = read_excel(file.path(base_dir, experiment, "elife-33105-supp5-v1.xlsx"),
                   sheet = 5) %>%
  # remove cells that did not pass the quality filter
  filter(pass_filter) %>%
  filter(pseudotime_state != 'NA') %>%
  dplyr::rename(label = pseudotime_state) %>%
  mutate(cell_type = 'p.falciparum') %>%
  column_to_rownames(var = 'sample_id')

# subset mat
mat02 = mat2[, c('id', rownames(meta2))]
mat02 %<>% column_to_rownames(var = 'id')

sc2 = CreateSeuratObject(mat02, meta = meta2)

#### Lastly do P.falciparum gametocytes
# load the count data
mat3 = read_excel(file.path(base_dir, experiment, "elife-33105-supp6-v1.xlsx"),
                  sheet = 3)

# load meta data for this experiment
meta3 = read_excel(file.path(base_dir, experiment, "elife-33105-supp5-v1.xlsx"),
                   sheet = 6) %>%
  # remove cells that did not pass the quality filter
  filter(pass_filter) %>%
  dplyr::rename(label = lasonder) %>%
  mutate(cell_type = 'p.falciparum.game') %>%
  column_to_rownames(var = 'sample_id')

# subset mat
mat03 = mat3[, c('id', rownames(meta3))]
mat03 %<>% column_to_rownames(var = 'id')

sc3 = CreateSeuratObject(mat03, meta = meta3)

# merge and save
sc_final = merge(sc2, sc3)
saveRDS(sc_final, file = file.path(output_dir, paste0(experiment, "_falciparum.rds")))
