# Aran et al. Nat Immunol 2019
# E-MTAB-7142

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Aran2018
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# Singler object to get cell types and post-QC cells
wget "https://github.com/dviraran/SingleR/blob/master/manuscript_figures/FiguresData/GSE111664.RData"

# Get raw count matrices from GEO
wget -O GSE111664_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE111664&format=file"
tar -xf GSE111664_RAW.tar
