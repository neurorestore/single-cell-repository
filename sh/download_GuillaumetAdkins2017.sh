# Guillaumet-Adkins et al., Genome Biol 2017: fresh vs. cryopreserved cell lines
# GSE102827

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=GuillaumetAdkins2017
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get matrices
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_Colon_Cryo.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_Colon_Fresh.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_HEK293_Cryo_Exp_1.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_HEK293_Cryo_Exp_2.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_HEK293_Cryo_Exp_3.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_HEK293_Fresh_Exp_1.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_HEK293_Fresh_Exp_2.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_HEK293_Fresh_Exp_3.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_HEK293_Nitrogen_Exp_2.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_K562_Cryo_Exp_1.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_K562_Cryo_Exp_2.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_K562_Fresh_Exp_1.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_K562_Fresh_Exp_2.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_MDCK_Cryo.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_MDCK_Fresh.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_NIH3T3_Cryo.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_NIH3T3_Fresh.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_PBMC_Cryo.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_PBMC_Fresh.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_PDOX_Cryo_Exp_1.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_PDOX_Cryo_Exp_2.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE85nnn/GSE85534/suppl/GSE85534_PDOX_Fresh_Exp_1.tsv.gz"
