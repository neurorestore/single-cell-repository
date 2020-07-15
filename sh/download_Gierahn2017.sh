# Gierahn et al., Nat Methods 2017: TB stimulated macrophages
# GSE92495

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Gierahn2017
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# raw data
wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92495/suppl/GSE92495_MTB.txt.gz'
