# McGinnis et al., Nat Methods 2019: 96-plex of HMECs, 15 different stimuli
# GSE129578

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=McGinnis2019
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129578/suppl/GSE129578_HMEC_orig_matrix.mtx.gz
## journal
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0433-8/MediaObjects/41592_2019_433_MOESM9_ESM.xlsx
