# Goldfarbmuren et al., Nat Commun 2020: lung epithelial cells after smoking
# GSE134174

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Goldfarbmuren2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get data
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134174/suppl/GSE134174_invivo_expression_matrix.txt.tar.gz'
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134174/suppl/GSE134174_Processed_invivo_metadata.txt.gz'

tar -xzf GSE134174_invivo_expression_matrix.txt.tar.gz
