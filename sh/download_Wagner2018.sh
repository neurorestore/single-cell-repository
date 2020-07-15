# Wagner et al., Science 2018: tyrosinase vs. chordin KO zebrafish
# GSE112294

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Wagner2018
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# raw data tarball
wget -O 'GSE112294_RAW.tar' 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE112294&format=file'
tar -xvf GSE112294_RAW.tar

# cluster names
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112294/suppl/GSE112294_ClusterNames.csv.gz
