## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: seurat_run_prostate.R
##
## Description: 
#~
## Wrapper for the seurat analysis of prostate data
##
## Authors: 
#~
## Stephan Gruener & Iros Barozzi
##
## License: 
#~
## GNU GPL v3
## Copyright 2024 
## Copyright Iros Barozzi Stephan Gruener
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Notes:
#~
##
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

options(bitmapType='cairo')

## Preparation

source("seurat_header.R")

## Technical parameters
n_cores <- 8                 # number of cores to use for parallel jobs
clustering_res <- 0.8        # resolution for clustering algorithm
sampling_N <- NA             # number of cells to randomly subsample, set NA if subsampling is not required

## Quality filtering thresholds (quite lenient, to accommodate all samples)
nFeature_RNA_min <- 250      # min number of genes expressed per cell
nCount_RNA_min <- 1000       # min number of transcripts expressed per cell
nCount_RNA_max <- 20000      # max number of transcripts expressed per cell
percent.mt_max <- 20         # max % of mito reads

## Metadata (full)
metadata_F <- "../data/metadata_muw_server.xlsx"
metadata_full <- read.xls(metadata_F, sheet = 1, header = TRUE, stringsAsFactors = FALSE)
# add id including sample & batch
metadata_full$id <- paste(metadata_full$description, metadata_full$batch, sep = "_batch-")

metadata_full_cmp <- read.xls(metadata_F, sheet = 2, header = TRUE, stringsAsFactors = FALSE)

#~~~~~~~~~~~~~~~~
#full data; LNCaP
#################

#~ LNCaP only

## Subsampling metadata, only run LNCaP and exclude WM-2d

# select LNCaP only
w <- grep("^LNCaP", metadata_full$description)
metadata <- metadata_full[w,]

w <- grep("^LNCaP", metadata_full_cmp$treated)
metadata_cmp <- metadata_full_cmp[w,]

# exclude 2d acute treatment
w <- grep("LNCaP_WM-2d", metadata$description, invert = T)
metadata <- metadata[w,]

w <- grep("LNCaP_WM-2d", metadata_cmp$treated, invert = T)
metadata_cmp <- metadata_cmp[w,]

## Using 5,000 features
n_features <- 5000

## Working Folder
outFolder <- paste("./", n_features/1000, "k_feat/", sep = "")
cmd <- paste("mkdir", outFolder)
system(cmd)

## Data preparation
source("seurat_data-prep.R")
print("Data-prep successful!")
gc() # clean environment

## Differential gene expression and cell cycle analysis
source("seurat_DEGs_cc.R")
print("Differential gene expression and cell cycle analysis successful!")
gc() # clean environment

source("seurat_SCTransform.R")
print("Normalization and module score calculations successful!")
#############################