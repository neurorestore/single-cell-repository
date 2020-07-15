# Mizoguchi et al., Nat Commun 2018
# GSE109449

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Mizoguchi2018
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109449/suppl/GSE109449_singlecell_rnaseq_gene_counts.tsv.gz"
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109449/suppl/GSE109449_singlecell_rnaseq_metadata.tsv.gz"
