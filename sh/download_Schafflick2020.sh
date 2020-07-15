# Schafflick et al., Nat Commun 2020: PBMCs and CSF in multiple sclerosis
# GSE138266

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Schafflick2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get UMI count matrices
wget -O GSE138266_RAW.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE138266&format=file'
tar -xvf GSE138266_RAW.tar

# get meta data
wget 'https://github.com/chenlingantelope/MSscRNAseq2019/blob/master/Notebooks/meta/celllabels.npy'
wget 'https://github.com/chenlingantelope/MSscRNAseq2019/blob/master/Notebooks/meta/doublet_score.npy'
wget 'https://github.com/chenlingantelope/MSscRNAseq2019/blob/master/Notebooks/meta/isMS.npy'
wget 'https://github.com/chenlingantelope/MSscRNAseq2019/blob/master/Notebooks/meta/predicted_doublets.npy'

# Note, converted cell type labels from pickle manually
