# Schirmer et al., Nature 2019: human brain in multiple sclerosis

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Schirmer2019
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# raw data
wget https://cells.ucsc.edu/ms/rawMatrix.zip
wget https://cells.ucsc.edu/ms/meta.tsv
unzip rawMatrix.zip
