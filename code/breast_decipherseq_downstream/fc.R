## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: fc.R
##
## Description: 
#~
## Investigate gene factor loadings 
## (as from iNMF/DECIPHER-seq)
## in different conditions
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


source("fc.inc.libs.R")

source("../lib/inc.func.R")

out_folder <- "fc_results/"
cmd <- paste0("mkdir ", out_folder)
system(cmd)

out_folder_traditiom <- "fc_results_traditiom/"
cmd <- paste0("mkdir ", out_folder_traditiom)
system(cmd)

# load metadata
metadata_F <- "../data/metadata_muw_server.xlsx"
metadata_full <- as_tibble(read.xls(metadata_F, sheet = 1, header = TRUE, stringsAsFactors = FALSE)) %>%
  mutate(sample = paste0(description, "_batch-", batch))

# load colors
colors <- read_tsv("../data/colors_for_umap.txt")

# load decipher-seq results
nmf_res <- readRDS("NMF_results_atK.rds")
rn <- tibble(rownames = rownames(nmf_res$MCF7$H))
name_tibble <- rn %>% 
  extract(rownames, into = c("sample", "barcode"), "(.*)_([^_]+)$", remove = FALSE) %>%
  left_join(metadata_full, by = c("sample" = "sample")) %>%
  dplyr::select(-stats_folder, -stats_file)

# join metadata with nmf results
H.data <- data.frame(nmf_res$MCF7$H) %>% 
  as_tibble(rownames = "rownames") %>%
  extract(rownames, into = c("sample", "barcode"), "(.*)_([^_]+)$", remove = FALSE) %>%
  left_join(metadata_full, by = c("sample" = "sample")) %>%
  left_join(colors, by = c("description" = "description")) %>%
  dplyr::select(-stats_folder, -stats_file) %>%
  relocate(starts_with("R"), .after = last_col())

# long format
H.data.long <- H.data %>%
  pivot_longer(cols = starts_with("R33"), names_to = "gene_program")

# prepare factors
source("fc.inc.prepare_factors.R")

# prepare gene sets for GSE/GSEA analyses
source("fc.inc.prepare_genesets.R")

# which parts of the script to run
run_clustering <- TRUE
run_umaps <- TRUE
run_factors <- TRUE
run_factors_sc <- TRUE
run_annotation <- TRUE
run_followup <- TRUE


#~~~~~~~~~~~
# Clustering
############

if (run_clustering) {

  # build SNN graph

  # Leiden clustering

  # assessing optimal resolution
  
}


#~~~~~~
# UMAPs
#######

if (run_umaps) {
  
  source("fc.run.umap.R")
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analyses of factor usage (global stats)
#########################################

if (run_factors) {
  
  source("fc.run.factor_usage.R")
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analyses of factor usage (single-cell)
########################################

if (run_factors_sc) {
  
  #note: sort the programs using the 'programs_order' variable
  #      (from the 'run_factors' step) to sort the results
  
  source("fc.run.factor_usage_sc.R")
  
}


#~~~~~~~~~~~~~~~~~~
# Factor annotation
###################

if (run_annotation) {
  
  #note: sort the programs using the 'programs_order' variable
  #      (from the 'run_factors' step) to sort the results
  
  source("fc.run.factor_anno.R")
  
  #TRADITIOM overlap
  
  source("fc.run.factor_anno.traditiom.R")
  
}


#~~~~~~~~~~~~~~~~~~
# Factor follow-ups
###################

if (run_followup) {
  
  source("fc.run.follow-up.IFN.R")
  
  source("fc.run.follow-up.survival.R")
  
}

