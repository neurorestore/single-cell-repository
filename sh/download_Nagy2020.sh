# Nagy et al., Nat Neurosci 2020: MDD vs. controls
# GSE141552

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Nagy2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get UMI count matrices
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz'

# get genes and cell type annotations
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_GeneNames.csv.gz'
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144136/suppl/GSE144136_CellNames.csv.gz'
