# Haber et al., Nature 2017: Salmonella treated endothelial cells
# GSE92332

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Haber2017
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# SalmHelm droplet
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_SalmHelm_UMIcounts.txt.gz"
