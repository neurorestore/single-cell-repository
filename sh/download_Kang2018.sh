# Kang et al., Nat Biotechnol 2018: IFN-b stimulated PBMCs
# GSE96583

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Kang2018
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# raw data
wget -O GSE96583_RAW.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE96583&format=file'
tar -xvf GSE96583_RAW.tar

# genes
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.genes.tsv.gz"

# meta
wget "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl/GSE96583_batch2.total.tsne.df.tsv.gz"
