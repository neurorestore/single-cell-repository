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
experiment = "Engelbrecht2020"

# read in loom data
loom = connect(file.path(base_dir, experiment,
                         "GSM4128644_scRNAseq.Aorta_GFP_EC.loom.data"))
mat = loom[['matrix']][,]
genes = loom[["row_attrs/Gene"]][]
cells = loom[["col_attrs/CellID"]][]
dimnames(mat) = list(cells, genes)
mat %<>% t()

# get meta data
meta = read_excel(file.path(base_dir, experiment,
  "GSM4128644_aorta.GFP.scRNAseq.cell.info.updated.xlsx")) %>%
  dplyr::rename(barcode = `Cell ID (Bam Filename)`) %>%
  mutate(barcode = paste0("Aorta_GFP_EC:", barcode)) %>%
  # drop missing cells
  drop_na() %>%
  dplyr::rename(label = Group)

# load in a manual download from their Pagoda bin
anno = readLines(file.path(base_dir, experiment, "GSE139065_download.txt"))
anno = seq(anno) %>%
  map(~ {
    split = unlist(strsplit(anno[.], ","))
    cell_type = split[2]
    cells = split[3:length(split)]
    out = data.frame(barcode = cells) %>%
      mutate(cell_type = cell_type)
    return(out)
  }) %>%
  bind_rows()

# bind this into meta data
meta %<>% left_join(anno)

# assign rownames
meta %<>% column_to_rownames('barcode')
mat %<>% extract(, rownames(meta))

# now back into Seurat
sc = CreateSeuratObject(mat, min.cells = 3, min.features = 0,
                        meta.data = meta)
# save
saveRDS(sc, file = file.path(output_dir, paste0(experiment, ".rds")))
