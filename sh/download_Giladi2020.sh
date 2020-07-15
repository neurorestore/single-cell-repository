# Giladi et al., Nat Immunol 2020: effect of EAE (MS) on monocytes
# GSE129150

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Giladi2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# count matrices
wget -O GSE144317_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144317&format=file"
tar -xvf GSE144317_RAW.tar

# clusters
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144317/suppl/GSE144317_sc_annotations.txt.gz

# sample level phenotype
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144317/suppl/GSE144317_metadata.txt.gz
