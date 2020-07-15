# Part 1: requires Seurat version 1.4 (and R version 3.4.3)
# script obtained from Supplementary Data of doi: 10.1038/s41586-018-0449-8
setwd("~/git/single-cell-repository")
options(stringsAsFactors = F)
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(tidyr)

# define base directory and output directory
base_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/raw")
output_dir = file.path(Sys.getenv("PROJECT"), "single-cell-repository/rnaseq/seurat")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

## write output files to a temporary folder
output_dir %<>% file.path('part_1')

# define experiment
experiment = "OrdovasMontanes2018"

# begin code from Ordovas-Montanes et al., Supplementary Data
polyp.data=read.table(file.path(base_dir, id, "SuppTable2_PolypALL_merged.txt"),
                      header=TRUE, row.names = 1)

#Initialize seurat object (names.delim set to '_' and Polyp1, Polyp2, etc. refers to biopsy taken, see Supplementary Table 9 for patient characteristics)
polyp <- new("seurat", raw.data = polyp.data)
polyp <- Setup(polyp, min.cells = 5, min.genes = 300, do.logNormalize = T,
               do.scale= T, do.center = T, names.field = 1, names.delim = "_",
               project = "Polyp_ALL_TOT")

#Add metadata to object for polyp status
polyp <- AddMetaData(polyp, NA, col.name = 'polyp')
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp1TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp2TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp3TOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp4TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp5TOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp6ATOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp6BTOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp7TOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp8TOT'] <- 'NO'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp9TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp11TOT'] <- 'YES'
polyp@data.info$polyp[polyp@data.info$orig.ident == 'Polyp12TOT'] <- 'YES'

#Filter cells with UMI counts greater than 12000
polyp <- SubsetData(polyp, subset.name = "nUMI", accept.high = 12000)

#Determine variable genes for input into PCA
polyp <- MeanVarPlot(polyp ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.13, x.high.cutoff = 7, y.cutoff = 0.28, do.contour = T)

#Run PCA
polyp <- PCA(polyp, pc.genes = polyp@var.genes, do.print = TRUE, pcs.print = 5,
             genes.print = 10)

#Find clusters and run TSNE over the first 12 PCs
polyp <- FindClusters(polyp, pc.use = 1:12, resolution = 1.2, print.output = 0,
                      save.SNN = T)

#Identify marker genes for clusters using a ROC test and print table
polyp.markers.roc= FindAllMarkers(polyp,return.thresh = 0.65, only.pos = TRUE,
                                  test.use = "roc", do.print = TRUE)

#Name cells for Figure 1
polyp <- AddMetaData(polyp, NA, col.name = 'subset')
polyp@data.info$subset[polyp@data.info$res.1.2 == 12] <- 'Basal'
polyp@data.info$subset[polyp@data.info$res.1.2 == 2] <- 'Basal'
polyp@data.info$subset[polyp@data.info$res.1.2 == 8] <- 'Basal'
polyp@data.info$subset[polyp@data.info$res.1.2 == 1] <- 'Apical'
polyp@data.info$subset[polyp@data.info$res.1.2 == 0] <- 'Apical'
polyp@data.info$subset[polyp@data.info$res.1.2 == 4] <- 'Apical'
polyp@data.info$subset[polyp@data.info$res.1.2 == 13] <- 'Glandular'
polyp@data.info$subset[polyp@data.info$res.1.2 == 3] <- 'Glandular'
polyp@data.info$subset[polyp@data.info$res.1.2 == 21] <- 'DOUBLET'
polyp@data.info$subset[polyp@data.info$res.1.2 == 16] <- 'Ciliated'
polyp@data.info$subset[polyp@data.info$res.1.2 == 9] <- 'TCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 19] <- 'DOUBLET'
polyp@data.info$subset[polyp@data.info$res.1.2 == 20] <- 'DOUBLET'
polyp@data.info$subset[polyp@data.info$res.1.2 == 15] <- 'PlasmaCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 17] <- 'PlasmaCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 7] <- 'PlasmaCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 10] <- 'PlasmaCell'
polyp@data.info$subset[polyp@data.info$res.1.2 == 5] <- 'Fibroblast'
polyp@data.info$subset[polyp@data.info$res.1.2 == 14] <- 'Fibroblast'
polyp@data.info$subset[polyp@data.info$res.1.2 == 6] <- 'Endothelial'
polyp@data.info$subset[polyp@data.info$res.1.2 == 11] <- 'Myeloid'
polyp@data.info$subset[polyp@data.info$res.1.2 == 18] <- 'MastCell'

#Removing biological cell doublets detected from analysis of marker gene lists and plot Figure 1B
select0 <- names(polyp@ident[polyp@ident == 0])
select1 <- names(polyp@ident[polyp@ident == 1])
select2 <- names(polyp@ident[polyp@ident == 2])
select3 <- names(polyp@ident[polyp@ident == 3])
select4 <- names(polyp@ident[polyp@ident == 4])
select5 <- names(polyp@ident[polyp@ident == 5])
select6 <- names(polyp@ident[polyp@ident == 6])
select7 <- names(polyp@ident[polyp@ident == 7])
select8 <- names(polyp@ident[polyp@ident == 8])
select9 <- names(polyp@ident[polyp@ident == 9])
select10 <- names(polyp@ident[polyp@ident == 10])
select11 <- names(polyp@ident[polyp@ident == 11])
select12 <- names(polyp@ident[polyp@ident == 12])
select13 <- names(polyp@ident[polyp@ident == 13])
select14 <- names(polyp@ident[polyp@ident == 14])
select15 <- names(polyp@ident[polyp@ident == 15])
select16 <- names(polyp@ident[polyp@ident == 16])
select17 <- names(polyp@ident[polyp@ident == 17])
select18 <- names(polyp@ident[polyp@ident == 18])
select19 <- names(polyp@ident[polyp@ident == 19])
select20 <- names(polyp@ident[polyp@ident == 20])
select21 <- names(polyp@ident[polyp@ident == 21])
nodoublet <- c(select0, select1, select2, select3, select4, select5, select6, select7, select8, select9, select10, select11, select12, select13, select14, select15, select16, select17, select18)
polyp <-SubsetData(polyp, cells.use = nodoublet)
polyp <- SetAllIdent(polyp,'subset')
polyp@ident <- factor(polyp@ident, levels(polyp@ident)[c(2,1,6,3,5,4,9,10,8,7)])

# write output
data = polyp@raw.data
meta = polyp@data.info
data = data %>%
  dplyr::select(rownames(meta))
saveRDS(data, file = file.path(output_dir, paste0(id, "_data.rds")))
saveRDS(meta, file = file.path(output_dir, paste0(id, "_meta.rds")))
