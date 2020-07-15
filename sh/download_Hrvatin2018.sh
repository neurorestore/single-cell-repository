# Hrvatin et al., Nat Neurosci 2018: mouse visual cortex after light exposure
# GSE102827

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Hrvatin2018
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# counts
wget -O 'GSE102827_RAW.tar' 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102827&format=file'
tar -xvf GSE102827_RAW.tar

# cell types
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102827/suppl/GSE102827_cell_type_assignments.csv.gz"
