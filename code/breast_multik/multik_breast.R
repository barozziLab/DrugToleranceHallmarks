## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: multik_breast.R
##
## Description: 
#~
## Optimal clustering resolution determination with multik
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

source("multiK_conserveMem.fct.R")
source("clustering.fct.R")

dataFolder <- "../breast_decipherseq_downstream/"
outFolder <- "./"
if(!dir.exists(outFolder)) dir.create(outFolder)

nmf_res <- readRDS(paste0(dataFolder, "NMF_results_atK.rds"))
H.matrix <- nmf_res$MCF7$H

rm(nmf_res)
gc()

maxMemory <- 500 # maximum memory to be allocated in GB
maxCores <- 120  # maximum number of cores to be used


reps = 100 
pSample = 0.1
cores = min(floor(maxMemory / (pSample * 50)), maxCores)
message(paste0("Running multiK with ", cores, " CPUs, ", reps, " reps and pSample of ", pSample, "."))

res <- MultiK_par_conserveMem(input.matr = H.matrix,
                              reps = reps,
                              pSample = pSample,
                              numCores_first = cores,
                              numCores_second = 120,
                              numChunks = 40)
saveRDS(res, file = paste0("results/grid_test_res_", pSample, ".rds"))

rm(res)
gc()

res <- readRDS(paste0("multik_res_breast_", pSample, ".rds"))

pdf(file = paste0("multik_res_breast_", pSample, ".pdf"), width = 18, height = 7)
res$plots
dev.off()






