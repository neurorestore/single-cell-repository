# Aztekin et al., Science 2019: limb amputation in regeneration-competent and incompetent tadpoles, and uninjuerd limbs
# E-MTAB-7716

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Aztekin2019
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7716/E-MTAB-7716.processed.1.zip

unzip E-MTAB-7716.processed.1.zi
unzip arrayExpressUpload.zip
unzip ArrayExpressV2.zip

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129578/suppl/GSE129578_POC_nuc_matrix.mtx.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129578/suppl/GSE129578_processed_data_files.tsv.tar.gz
tar -xzvf GSE129578_processed_data_files.tsv.tar.gz

## journal
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-019-0433-8/MediaObjects/41592_2019_433_MOESM3_ESM.xlsx
