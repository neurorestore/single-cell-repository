# Engelbrecht et al., eLife 2020; Sphingosine 1-phosphate-regulated transcriptomes

## define your PROJECT directory
PROJECT=${PROJECTS}/single-cell-repository

EXPERIMENT=Engelbrecht2020
mkdir -p ${PROJECT}/rnaseq/raw/${EXPERIMENT}
cd ${PROJECT}/rnaseq/raw/${EXPERIMENT}

# get matrices and meta data
wget -O GSE139065_RAW.tar 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139065&format=file'
tar -xf GSE139065_RAW.tar

# bulk data
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139065/suppl/GSE139065_Supplementary.File.1.Bulk.RNA-seq.GFP.and.ECKO.xls.gz'
gunzip GSE139065_Supplementary.File.1.Bulk.RNA-seq.GFP.and.ECKO.xls.gz

# To get the cell type annotations:
# Download their Binary ‘.bin’ file that can be uploaded to the graphical user interface: at:
# https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNTI2OTAvZWxpZmUtNTI2OTAtc3VwcDctdjIuYW9ydGEuYmlu/elife-52690-supp7-v2.aorta.bin?_hash=4cpUWf6x7%2FODC203ChKnW6G0LNuiJCYzz59suWpqwVA%3D

# Download the annoations at:
# https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvNTI2OTAvZWxpZmUtNTI2OTAtc3VwcDgtdjIuYW9ydGE=/elife-52690-supp8-v2.aorta?_hash=it%2BCGMOM2QfDo4h1fGnq6QO1E%2FdkK%2FDmrXRpget9TeU%3D

# Load this into the http://pklab.med.harvard.edu/nikolas/pagoda2/frontend/current/pagodaLocal/.
# Export the relevant cell types
# Load this back into R
# Parse the  output and bind to other provided meta data
