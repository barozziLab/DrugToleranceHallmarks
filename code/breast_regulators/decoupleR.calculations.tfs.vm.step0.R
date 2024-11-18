
## Notes
## Runs on the MUW virtual machine
## First load the right conda environment: 
## --> conda activate decoupler_202407

## Folder:
## --> cd /home/iros/HDBC/decoupleR_full/

library(tidyverse)
library(ggplot2)
library(Seurat)

so <- readRDS("so_cca_integrated_cc.rds")

gene_table <- read_tsv("features.tsv.gz", col_names = c("ENSID", "Symbol", "Anno"))

cell_types <- so@meta.data$cell_type %>% unique() %>% sort()

for (cell_type_oi in cell_types) {

  so_subset <- subset(x = so, subset = cell_type == cell_type_oi)

  outF <- paste0("/share/iros/input.", cell_type_oi, ".rds")

  # Save res object for downstream analyses
  saveRDS(object = so_subset, file = outF)

  rm(so_subset)
  gc()

}
