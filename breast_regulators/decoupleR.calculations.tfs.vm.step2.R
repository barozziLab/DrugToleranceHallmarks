
## Notes
## Runs on the MUW virtual machine
## First load the right conda environment: 
## --> conda activate decoupler_202407

## Folder:
## --> cd /home/iros/HDBC/decoupleR_full/

library(tidyverse)
library(Seurat)
library(decoupleR)


#~~~~~~~~~~~~~~~~~~~~~~~
#Summarise TF-activities
########################

pvalue_thresh <- 0.05

cell_types <- gsub("\\.rds", "", gsub("input\\.", "", system("ls /share/iros", intern = TRUE)))

tf_stats <- list()

for (cell_type_oi in cell_types) {

    print(paste0("Processing ", cell_type_oi, " ... "))

    inF <- paste0("decoupleR.activities.tfs.split.", cell_type_oi, ".rds")
    net_tf_acts <- readRDS(inF)

    #TFs to be considered (= at least one significant association)
    tfs_to_consider <- net_tf_acts %>% 
        dplyr::filter(p_value <= pvalue_thresh) %>% 
        pull(source) %>% 
        unique()

    #Statistics per TF (to be considered, see above)
    tf_stats[[cell_type_oi]] <- net_tf_acts %>% 
        mutate(significant = p_value <= pvalue_thresh) %>%
        dplyr::filter(source %in% tfs_to_consider) %>% 
        mutate(cell_type = cell_type_oi) %>%
        group_by(source, cell_type) %>%
        summarise(mean = mean(score), n_cells_pos = sum(significant), n_cells_tot = n()) %>%
        mutate(n_cell_frac = n_cells_pos / n_cells_tot) %>% 
        ungroup()

    rm(net_tf_acts)
    gc()

}

#merge
tf_stats_merge <- do.call(rbind, tf_stats)

#save
saveRDS(object = tf_stats_merge, file = "decoupleR.activities.tfs.merge.tf_stats_merge.rds")


#~~~~~~~~~~~~~~~~~~
#TF-gene expression
###################

so <- readRDS("so_cca_integrated_cc.rds")

# scTransform normalized data (integrated)
mat <- as.matrix(so@assays$integrated@data)

tf_ids <- sort(unique(tf_stats_merge$source))

df_tf_gene <- t(mat[rownames(mat) %in% tf_ids,]) %>%
  as.data.frame() %>%
  mutate(cell_type = so@meta.data$cell_type) %>%
  pivot_longer(cols = -c(cell_type), names_to = "TF_gene", values_to = "expression")

# Fraction of cells expressing the gene per cell type
tf_expr_frac <- df_tf_gene %>%
  group_by(cell_type, TF_gene) %>%
  summarise(more0 = sum(expression > 0), tot = n()) %>%
  mutate(fraction = more0 / tot)

#save
saveRDS(object = tf_expr_frac, file = "decoupleR.activities.tfs.merge.tf_expr_frac.rds")

