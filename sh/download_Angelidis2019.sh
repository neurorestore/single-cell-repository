# Angelidis et al., Nat Commun 2019: aging mouse lung
# GSE124872

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Angelidis2019
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124872/suppl/GSE124872_raw_counts_single_cell.RData.gz'
gunzip GSE124872_raw_counts_single_cell.RData.gz

## follow the steps from : https://github.com/theislab/2018_Angelidis/blob/master/download.sh
wget https://hmgubox.helmholtz-muenchen.de/f/c998398ccff04dd39237/?dl=1 -O 2018_Angelidis.tar
mkdir data; cd data
tar -xvf ../2018_Angelidis.tar
