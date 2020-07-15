# Hu et al., Mol Cell 2017: cortical tissues +/- pentylenetetrazole
# GSE106678

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Hu2017
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# raw data
# note: Methods refers to GitHub: https://github.com/wulabupenn/Hu_MolCell_2017
# GitHub points to Google Drive: https://drive.google.com/drive/folders/1g7g_eoXX3oA6i248Tdx53xGqz2ZF1Ga6
# download data.tar.gz from here

# untar
tar -xvf data.tar.gz
