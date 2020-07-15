# Farbehi et al., eLife 2019: mouse cardiac interstitial cells 3 & 7 days after sham or myocardial infarction injury
# E-MTAB-7376

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Farbehi2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

wget 'https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7376/E-MTAB-7376.processed.1.zip'
wget 'https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7376/E-MTAB-7376.processed.2.zip'
wget 'https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7376/E-MTAB-7376.processed.3.zip'

unzip E-MTAB-7376.processed.1.zip
unzip E-MTAB-7376.processed.2.zip
unzip E-MTAB-7376.processed.3.zip
