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
experiment = "Hagai2018"

# meta
files = list.files(file.path(base_dir, experiment), pattern = "*.txt.gz",
  full.names = T)

# do each species independently due to different gene usage
species_list = c("rat", "rabbit", "pig", "mouse")

# define GTF directory to map ENSEMBL to SYMBOL
gtf_dir = "~/projects/rrg-aphil/aphil/sclq/ensembl"

for (species in species_list) {
  # get GTF and map ENSEMBL to SYMBOL
  if (species == 'mouse') {
    gtf = 'Mus_musculus.GRCm38.95.gtf'
  } else if (species == 'rat') {
    gtf = 'Rattus_norvegicus.Rnor_6.0.95.gtf'
  } else if (species == 'rabbit') {
    gtf = 'Oryctolagus_cuniculus.OryCun2.0.99.gtf'
  } else if (species == 'pig') {
    gtf = 'Sus_scrofa.Sscrofa11.1.99.gtf'
  }
  # get GTF
  gtf = as.data.frame(rtracklayer::import(file.path(gtf_dir, gtf)))

  sub_files = files[grepl(species, files)]
  # Take the first matrix as a ref for gene order
  gene_ref = fread(sub_files[1]) %>% pull(V1)
  # initialise final matrix
  mat2 = data.frame(gene = gene_ref)
  meta2 = data.frame()
  for (i in seq(1, length(sub_files))){
    file = sub_files[i]
    # Load in data
    mat = fread(file) %>%
      column_to_rownames(var="V1")
    col_names = paste0(i, "-", colnames(mat))
    mat %<>% set_colnames(col_names)
    # Add to final matrix
    mat2 = cbind(mat2, mat)
    # get meta data
    meta = data.frame(
      sample = rep(basename(file), ncol(mat))
    ) %>%
      separate(sample, c("replicate", "label", "filtered", "by",
        "cell", "cluster"), "_", remove = F) %>%
      select(-filtered, -by, -cell, -cluster) %>%
      mutate(cell_type = 'bone marrow derived mononuclear phagocytes') %>%
      set_rownames(colnames(mat))
    meta2 = rbind(meta2, meta)
  }
  mat2 %<>% select(-gene)

  ensembl_map = gtf %>%
    dplyr::select(gene_id, gene_name) %>%
    dplyr::rename(ensembl = 'gene_id', symbol = 'gene_name') %>%
    dplyr::distinct() %>%
    filter(!is.na(symbol)) %>%
    filter(ensembl %in% rownames(mat2)) %>%
    # take the first match
    group_by(symbol) %>%
    dplyr::slice(1)

  # subset expr and set colnames
  mat2 %<>% extract(ensembl_map$ensembl,) %>% set_rownames(ensembl_map$symbol)

  # create Seurat object
  sc = CreateSeuratObject(mat2, min.cells = 3, min.features = 0,
    meta.data = meta2)
  # save
  saveRDS(sc, file = file.path(output_dir, paste0(experiment, "_", species, ".rds")))
}
