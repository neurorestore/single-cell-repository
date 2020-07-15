# Wang et al., Cell 2020: aging primate ovary
# GSE130664

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Wang2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get UMI count matrices barcode map
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130664/suppl/GSE130664_barcode_information.txt.gz'
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130664/suppl/GSE130664_merge_UMI_count.txt.gz'

# cell metadata obtained from Table S1 of manuscript
