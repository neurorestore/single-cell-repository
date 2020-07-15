# Hagai et al., Nature 2018: bone marrow derived mononuclear phagocytes derived from mouse; rat; rabbit or pig

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Hagai2018
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# Second zip directory which contains UMI matrices
wget "https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6754/E-MTAB-6754.processed.2.zip"
unzip E-MTAB-6754.processed.2.zip
