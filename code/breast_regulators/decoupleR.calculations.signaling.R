
## Notes
## Runs on the MUW virtual machine
## First load the right conda environment: 
## --> conda activate decoupler_202407

## Folder:
## --> cd /home/iros/HDBC/decoupleR_full/

library(tidyverse)
library(Seurat)
library(decoupleR)

so <- readRDS("so_cca_integrated_cc.rds")


#~~~~~~~~~~~~~
#Prepare Input
##############

# scTransform normalized data (integrated)
mat <- as.matrix(so@assays$integrated@data)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Compute Signalling Activities
##############################

net_sig <- get_progeny(organism = 'human', top = 500)

# Run mlm
net_sig_acts <- run_mlm(mat=mat, 
                        net=net_sig, 
                        .source='source', 
                        .target='target', 
                        .mor='weight', 
                        minsize = 5)

# Extract mlm and store it in a new assay called pathwaysmlm in the Seurat object
so[['pathwaysmlm']] <- net_sig_acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)


#~~~~~~~~~~~
#Save Output
############

acts_sig_path <- "decoupleR.activities.sig.rds"

#so_decoupleR_path <- "/share/iros/so_obj_cca_integrated_subset.20240723.rds"

# Activities
saveRDS(object = net_sig_acts, file = acts_sig_path)

# Save the updated object
#saveRDS(object = so, file = so_decoupleR_path)


