# Huang et al., Cell 2020: IBD in pediatric cohort

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Huang2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get files
wget "https://zhanglaboratory.com/wp-content/uploads/2020/02/Pediatric_IBD_colitis_scRNA_countMatrix.txt.zip"
wget "https://zhanglaboratory.com/wp-content/uploads/2020/02/Pediatric_IBD_colitis_scRNA_annotated.txt.zip"

unzip Pediatric_IBD_colitis_scRNA_annotated.txt.zip
unzip Pediatric_IBD_colitis_scRNA_countMatrix.txt.zip
