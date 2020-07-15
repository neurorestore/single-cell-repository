# Jakel et al., Nature 2019: human oligodendrocytes in multiple sclerosis
# GSE118257

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Jakel2019
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# raw data
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118257/suppl/GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE118nnn/GSE118257/suppl/GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz
