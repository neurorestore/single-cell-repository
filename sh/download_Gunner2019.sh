# Gunner et al., Nat Neurosci 2019: effect of whisker lesioning on the somatosensory cortex (synapse elimination) in WT or CX3CL1 KO mice
# GSE129150

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Gunner2019
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# count matrices
wget -O 'GSE129150_RAW.tar' https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE129150&format=file
tar -xzvf GSE129150_RAW.tar

# clusters
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129150/suppl/GSE129150_seuratClusters.csv.gz

# sample level phenotype
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129150/suppl/GSE129150.xls.gz
