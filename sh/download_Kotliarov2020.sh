# Kotliarov et al., Nat Med 2020: vaccine response in PBMCs

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Kotliarov2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get data
wget -o H1_day0_scranNorm_adtbatchNorm_dist_clustered_TSNE_labels.rds 'https://nih.figshare.com/ndownloader/files/20706645'
wget -o clustree_node_labels_withCellTypeLabels.txt 'https://nih.figshare.com/ndownloader/files/20706633'
