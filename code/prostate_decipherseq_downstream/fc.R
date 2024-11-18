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


source("../breast_decipherseq_downstream/fc.inc.libs.R")

source("../lib/inc.func.R")

source("../lib/inc.func.decoupleR.R")

out_folder <- "fc_results/"
cmd <- paste0("mkdir ", out_folder)
system(cmd)

# load metadata
metadata_F <- "../data/prostate/metadata_filtered_ordered.xlsx"
metadata_full <- as_tibble(read.xls(metadata_F, sheet = 2, header = TRUE, stringsAsFactors = FALSE)) %>%
  mutate(sample = paste0(description, "_batch-", batch))
# exclude LNCaP_WM-2d
metadata_full <- metadata_full %>% dplyr::filter(description != "LNCaP_WM-2d")

# load decipher-seq results
nmf_res <- readRDS("NMF_results_atK.rds")
rn <- tibble(rownames = rownames(nmf_res$LNCaP$H))
name_tibble <- rn %>% 
  extract(rownames, into = c("sample", "barcode"), "(.*)_([^_]+)$", remove = FALSE) %>%
  left_join(metadata_full, by = c("sample" = "sample")) %>%
  select(-stats_folder, -stats_file)

# join metadata with nmf results
H.data <- data.frame(nmf_res$LNCaP$H) %>% 
  as_tibble(rownames = "rownames") %>%
  extract(rownames, into = c("sample", "barcode"), "(.*)_([^_]+)$", remove = FALSE) %>%
  left_join(metadata_full, by = c("sample" = "sample")) %>%
  select(-stats_folder, -stats_file) %>%
  relocate(starts_with("R"), .after = last_col())

# long format
H.data.long <- H.data %>%
  pivot_longer(cols = starts_with("R22"), names_to = "gene_program")

# which parts of the script to run
run_factors <- TRUE
run_annotation <- TRUE
run_primary <- TRUE
run_decoupleR <- TRUE

#~~~~~~~~~~~~~~~~~~~~~~~~~
# Analyses of factor usage
##########################

if (run_factors) {
  
  source("fc.run.factor_usage.R")
  
}

#~~~~~~~~~~~~~~~~~~
# Factor annotation
###################

if (run_annotation) {
  
  #note: sort the programs using the 'programs_order' variable
  #      (from the 'run_factors' step) to sort the results
  
  source("fc.run.factor_anno.R")
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Integration with primary data
###############################

if (run_primary) {
  
  source("fc.run.primary_data.R")
  
}

#~~~~~~~~~~~~~~~~~~~~~~~
# Regulators (decoupleR)
########################

if (run_decoupleR) {
  
  source("fc.run.decoupleR.main.R")
  
}

