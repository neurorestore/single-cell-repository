# Ordovas-Montanes et al., Nature 2018: patients with chronic rhinosinusitis with or without nasal polyps

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=OrdovasMontanes2018
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# raw data is in a supplementary table
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0449-8/MediaObjects/41586_2018_449_MOESM4_ESM.zip
unzip 41586_2018_449_MOESM4_ESM.zip
unzip SuppTable2_PolypALL_merged.txt.zip
