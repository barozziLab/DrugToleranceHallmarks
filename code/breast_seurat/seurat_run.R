## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: seurat_run.R
##
## Description: 
#~
## Wrapper for the seurat analysis of breast data
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
## Run on HPC
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

options(bitmapType='cairo')

## Preparation
#vsc_data_dir <- "/gpfs/data/fs72016/stephangrun/ET_resistance/seurat_breast_final_06_2024/"
working_dir <- "./"

source("seurat_header.R")

## Technical parameters
n_cores <- 8              # number of cores to use for parallel jobs
clustering_res <- 0.8     # resolution for clustering algorithm
sampling_N <- NA          # number of cells to randomly subsample, set NA if subsampling is not required

## Quality filtering thresholds (quite lenient, to accommodate all samples)
nFeature_RNA_min <- 250   # min number of genes expressed per cell
nCount_RNA_min <- 1000    # min number of transcripts expressed per cell
nCount_RNA_max <- 20000   # max number of transcripts expressed per cell
percent.mt_max <- 20      # max % of mito reads

## Metadata (full)
metadata_full <- as.data.frame(read_tsv("metadata_full.tsv"))
# add id including sample & batch
metadata_full$id <- paste(metadata_full$description, metadata_full$batch, sep = "_batch-")

metadata_full_cmp <- as.data.frame(read_tsv("metadata_full_cmp.tsv"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# full data; CCA kit
#############################

#~ MCF7 only
  
## Subsampling metadata, only run MCF7 (no LNCaP)
w <- setdiff(grep("^MCF7", metadata_full$description), grep("CD44", metadata_full$description))
metadata <- metadata_full[w,]

w <- setdiff(grep("^MCF7", metadata_full_cmp$treated), grep("CD44", metadata_full_cmp$treated))
metadata_cmp <- metadata_full_cmp[w,]

## Using 5,000 features
n_features <- 5000

## Working Folder
outFolder <- paste(working_dir, "/MCF7_cca_kit_", n_features/1000, "k_feat/", sep = "")
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

## CCA integration and 
source("seurat_cca-pairs.kit.R")
print("Integration by kit version successful!")
#############################