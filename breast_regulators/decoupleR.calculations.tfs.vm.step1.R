
## Notes
## Runs on the MUW virtual machine
## First load the right conda environment: 
## --> conda activate decoupler_202407

## Folder:
## --> cd /home/iros/HDBC/decoupleR_full/

library(tidyverse)
library(Seurat)
library(decoupleR)


#~~~~~~~~~~~~~~~~~~~~~
#Compute TF-activities
######################

net_tf <- get_collectri(organism='human', split_complexes=FALSE)

cell_types <- gsub("\\.rds", "", gsub("input\\.", "", system("ls /share/iros", intern = TRUE)))

for (cell_type_oi in cell_types) {

    print(paste0("Processing ", cell_type_oi, " ... "))

    inF <- paste0("/share/iros/input.", cell_type_oi, ".rds")
    so_subset <- readRDS(inF)

    mat <- as.matrix(so_subset@assays$integrated@data)
    rm(so_subset)
    gc()

    # Run ulm
    net_tf_acts <- run_ulm(mat=mat, 
                           net=net_tf, 
                           .source='source', 
                           .target='target',
                           .mor='mor', 
                           minsize = 10)

    outF <- paste0("decoupleR.activities.tfs.split.", cell_type_oi, ".rds")

    # Save res object for downstream analyses

    saveRDS(object = net_tf_acts, file = outF)

    rm(net_tf_acts)
    gc()

}

