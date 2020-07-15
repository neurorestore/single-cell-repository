setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(tidyverse)
library(Matrix)
library(Seurat)
library(magrittr)
library(data.table)
library(readxl)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# define experiment
experiment = "Gunner2019"

# read count matrices
files = list.files(file.path(base_dir, experiment),
                   full.names = T,
                   pattern = "GSM37*")
dats = map(files, fread)

# read clusters
clust = read.csv(file.path(base_dir, experiment,
                           "GSE129150_seuratClusters.csv.gz")) %>%
  set_colnames(c("barcode", "cell_type"))

# manually recode clusters based on figures in the paper
clust %<>%
  mutate(cell_type = fct_recode(as.character(cell_type),
                                '0_Excitatory - layer IV' = '0',
                                '1_Excitatory' = '1',
                                '2_Excitatory' = '2',
                                '3_Excitatory' = '3',
                                '5_Excitatory' = '5',
                                '13_Excitatory' = '13',
                                '15_Excitatory' = '15',
                                '17_Excitatory' = '17',
                                '19_Excitatory' = '19',
                                '23_Excitatory' = '23',
                                '7_Inhibitory' = '7',
                                '8_Inhibitory' = '8',
                                '9_Inhibitory' = '9',
                                '10_Inhibitory' = '10',
                                '11_Inhibitory' = '11',
                                '14_Inhibitory' = '14',
                                '18_Inhibitory' = '18',
                                '21_Inhibitory' = '21',
                                '16_Microglia' = '16',
                                '4_Astrocytes' = '4',
                                '6_Astrocytes' = '6',
                                '12_OPCs' = '12',
                                '22_OPCs' = '22',
                                '25_Oligodendrocytes' = '25',
                                '26_Macrophages' = '26',
                                '20_Endothelial cells' = '20',
                                '24_Pericytes' = '24',
                                '27_Unidentified' = '27'
                                ))

# read phenotypes
system(paste("gunzip", file.path(base_dir, experiment, "GSE129150.xls.gz")))
pheno = read_excel(file.path(base_dir, experiment, "GSE129150.xls"))

# merge datasets
samples = gsub("\\.counts.*$", "", basename(files)) %>%
  map(~ unlist(strsplit(., "_"))) %>%
  map(~ .[-1]) %>%
  map_chr(~ paste0(., collapse = "_"))
dat = dats %>%
  setNames(samples) %>%
  rbindlist(idcol = 'sample')

# merge in cluster labels
dat %<>%
  mutate(barcode = paste0(sample, "_", barcode)) %>%
  left_join(clust, by = 'barcode')

# merge in phenotypes
pheno %<>%
  dplyr::select(title,
                `characteristics: genotype`,
                `characteristics: sex`,
                `characteristics: treatment`) %>%
  mutate(title = gsub("Het_", "Het", title),
         title = gsub("KO_", "KO", title)) %>%
  set_colnames(c("sample", "genotype", "sex", "label"))
dat %<>% left_join(pheno, by = 'sample')

# extract combined metadata
meta = dat %>%
  dplyr::select(barcode, cell_type, sample, genotype, sex, label) %>%
  rename(replicate = "sample") %>%
  column_to_rownames('barcode')

# drop cells that were not included in Supp. Fig. 10
meta %<>% drop_na(cell_type)

# convert to gene-by-cell matrix
mat = dat %>%
  dplyr::select(-cell_type, -sample, -genotype, -sex, -label) %>%
  column_to_rownames('barcode') %>%
  as.matrix()
keep = rownames(meta)
mat %<>% extract(keep, )

# create Seurat object
sc = CreateSeuratObject(t(mat), min.cells = 3, min.features = 0,
  meta.data = meta)

# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
