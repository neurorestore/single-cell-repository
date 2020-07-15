# Mathys et al., Nature 2019: Alzheimer's brain tissue
# download from synapse

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Mathys2019
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

## install synapse python tools

# data
synapse get -r syn18681734

# meta data
synapse get -r syn18642926

# patient data
synapse get syn3191087

# meta data from paper
wget "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1195-2/MediaObjects/41586_2019_1195_MOESM5_ESM.xlsx"
