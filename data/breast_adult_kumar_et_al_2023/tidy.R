
library(tidyverse)

m2h <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/mouse.txt", col_names = c("gene", "ortholog"))
h2m <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/human.txt", col_names = c("gene", "ortholog"))

gsets <- c()

data_clean <- read_tsv("raw/41586_2023_6252_MOESM1_ESM.signatures.good_formatting.txt")

#add consolidated signature name
data_clean <- data_clean %>% 
  rowwise() %>%
  mutate(prefix = ifelse(table %in% c("S3_cell_types", "S4_cell_type_nuclei"), "cell-type_", "cell-state_")) %>%
  mutate(gs_name = paste0(prefix, signature)) %>%
  select(gs_name, gene)

#add mouse symbols
data_clean <- data_clean %>%
  left_join(h2m, by = "gene") %>%
  dplyr::filter(!is.na(ortholog))

#split and re-pool
data_clean_hs <- data_clean %>%
  select(gs_name, gene) %>%
  mutate(species = "hs") %>%
  unique()
data_clean_mm <- data_clean %>%
  select(gs_name, ortholog) %>%
  dplyr::rename(gene = ortholog) %>%
  mutate(species = "mm") %>%
  unique()
data_clean_pool <- rbind(data_clean_hs, data_clean_mm)

write_tsv(x = data_clean_pool, file = "markers.txt")

