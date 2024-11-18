## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: multik_prostate.R
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

dataFolder <- "../prostate_decipherseq_downstream/"
outFolder <- "./"
if(!dir.exists(outFolder)) dir.create(outFolder)

nmf_res <- readRDS(paste0(dataFolder, "NMF_results_atK_prostate_no_2d_acute.rds"))
H.matrix <- nmf_res$LNCaP$H

rm(nmf_res)
gc()

maxMemory <- 500 # maximum memory to be allocated in GB
maxCores <- 120  # maximum number of cores to be used
  
reps = 100 
pSample = 0.8
cores = 120

message(paste0("Running multiK with ", cores, " CPUs, ", reps, " reps and pSample of ", pSample, "."))

res <- MultiK_par_conserveMem(input.matr = H.matrix,
                              reps = reps,
                              pSample = pSample,
                              numCores_first = cores,
                              numCores_second = cores,
                              numChunks = 30)
saveRDS(res, file = paste0(outFolder, "multik_res_prostate_no_2d_acute_", pSample, ".rds"))
  
rm(res)
gc()

res <- readRDS(paste0(outFolder, "multik_res_prostate_no_2d_acute_", pSample, ".rds"))

pdf(file = paste0(outFolder, "multik_res_prostate_no_2d_acute_", pSample, ".pdf"), width = 18, height = 7)
res$plots
dev.off()





