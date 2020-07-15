# Wu et al., Neuron 2017: scRNA-seq with actinomycin D dissociation
# GSE103976

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Wu2017
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# raw data
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103976/suppl/GSE103976_MeA_AllCells_DGE.txt.gz

# neuronal DGE (for subtypes)
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103976/suppl/GSE103976_MeA_Neurons_DGE.txt.gz
