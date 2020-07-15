# Smillie et al., Cell 2019: ulcerative colitis
# data is stored on Broad Institute 'Single Cell Portal'

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Smillie2019
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# curl command on URLs from:
# "https://portals.broadinstitute.org/single_cell/study/SCP259/intra-and-inter-cellular-rewiring-of-the-human-colon-during-ulcerative-colitis#/"
