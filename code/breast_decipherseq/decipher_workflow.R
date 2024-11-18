## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: decipher_workflow.R
##
## Description: 
#~
## Run adapted DECIPHER-seq workflow on breast data
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
## Long runtime, high memory allocation. Run on VM/HPC
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# DECIPHER-seq workflow for one celltype
library(igraph)
library(Seurat)
library(ggplot2)
library(rliger)
library(RANN)
library(rlist)
library(ape)
library(matrixStats)
library(patchwork)
library(geiger)
library(tidyverse)

source("DECIPHER-seq_adapted/DECIPHER-seq.functions_v2.R")
source("DECIPHER-seq_adapted/DECIPHER-seq.util.R")

so <- readRDS(file = "../breast_seurat/MCF7_cca_kit_5k_feat/so.rds")

#~~~~~~~~~~
so[['Sample']] <- "Sample_MCF7"
so[['Type']] <- "MCF7"
so[['Batch']] <- so[['kit']]

Idents(so) <- 'kit'

# This step was only included for testing with a subset of cells
#so <- so[,WhichCells(so, idents = c("v2", "v3"), downsample = 500)]
#~~~~~~~~~~

# The results for this step can be loaded directly
# NMF_results <- iNMF_ksweep(so, k.max = 40, max.cores = 1)
# saveRDS(NMF_results, file = "results/NMF_results.rds")
NMF_results <- readRDS("NMF_results.rds")

phylo_trees <- lapply(NMF_results, build_phylo_tree)

program_outlier_score <- lapply(NMF_results, identify_outlier_programs)

suggested_thresholds = suggest_dist_thresh(phylo_trees)
thresh_use = round(max(suggested_thresholds), 1)

thresh_grid <- seq(0.1, 2, 0.1)

p = plot_hists(phylo_trees, thresh.use = thresh_use)

suggested_thresholds = suggest_dist_thresh(phylo_trees)
thresh_use = 0.3
phylo_partitions = mapply(partition_phylo_tree, x = phylo_trees, y = program_outlier_score, dist.thresh = thresh_use, outlier.thresh = 5, SIMPLIFY = F)

par(mfrow=c(1,5), mar=c(2.5,5.5,1,1))
plot_phylo_trees(phylo_trees, phylo_partitions)
  
K_metrics = lapply(phylo_partitions, calculate_K_metric)
k.use = lapply(K_metrics, suggest.k)
  
lineplot <- ggplot(K_metrics[[1]] %>% mutate(rank_n = as.integer(str_sub(rank, start = 2))), aes(x = rank_n, y = weighted_n_subtrees)) +
    geom_path(color = "darkgreen") +
    geom_point(color = "darkgreen") +
    geom_vline(xintercept = k.use[[1]]) +
    ggtitle(paste0("threshold: ", thresh_use)) +
    theme_bw()

write_tsv(lineplot$data, file = "k_metrics_lineplot.tsv")

ggsave("lineplot.pdf", lineplot, width = 6, height = 2)

NMF_results_atK <- mapply(FUN = NMF_results_opt_k, NMF_results, k.use, program_outlier_score, SIMPLIFY = F)
saveRDS(NMF_results_atK, file = "../breast_decipherseq_downstream/NMF_results_atK.rds")
