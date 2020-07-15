# Bhattacherjee et al., Nat Commun 2019: PFC of mice withdrawing from cocaine
# GSE124952

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Bhattacherjee2019
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get files
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124952/suppl/GSE124952_expression_matrix.csv.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124952/suppl/GSE124952_meta_data.csv.gz"
