# Reid et al., eLife 2020: scRNA-seq of malaria parasites

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Reid2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get meta data
wget -O elife-33105-supp5-v1.xlsx 'https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzMxMDUvZWxpZmUtMzMxMDUtc3VwcDUtdjEueGxzeA==/elife-33105-supp5-v1.xlsx?_hash=rTNIWPO0PuW8yvwWPc%2BYNOgmqugLWT6wumJQJnmZDIA%3D'

# get matrices
wget -O elife-33105-supp6-v1.xlsx 'https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzMxMDUvZWxpZmUtMzMxMDUtc3VwcDYtdjEueGxzeA==/elife-33105-supp6-v1.xlsx?_hash=qFCrRgVcefi0KJsUMW07In65jzVe0TUCQedjJk4CuY4%3D'
