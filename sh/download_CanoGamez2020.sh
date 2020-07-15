# Cano-Gamez et al., Nat Commun 2020: scRNA-seq with cytokine induced T cell differentiation

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=CanoGamez2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get matrix information
wget 'https://storage.googleapis.com/otar040-cytoimmunogenomics/publication-project-data/scRNAseq/NCOMMS-19-7936188_scRNAseq_barcodes.tsv'
wget 'https://storage.googleapis.com/otar040-cytoimmunogenomics/publication-project-data/scRNAseq/NCOMMS-19-7936188_scRNAseq_genes.tsv'
wget 'https://storage.googleapis.com/otar040-cytoimmunogenomics/publication-project-data/scRNAseq/NCOMMS-19-7936188_scRNAseq_raw_UMIs.mtx'

# get meta data
wget 'https://storage.googleapis.com/otar040-cytoimmunogenomics/publication-project-data/scRNAseq/NCOMMS-19-7936188_metadata.txt'
