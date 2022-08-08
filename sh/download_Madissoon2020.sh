# Madissoon et al., Genome Biol 2020: scRNA-seq with cold preservation in 3 tissues

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Madissoon2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get data from https://www.tissuestabilitycellatlas.org/
wget 'https://cellgeni.cog.sanger.ac.uk/tissue-stability/oesophagus_ts.rds'
wget 'https://cellgeni.cog.sanger.ac.uk/tissue-stability/lung_ts.rds'
wget 'https://cellgeni.cog.sanger.ac.uk/tissue-stability/spleen_ts.rds'
