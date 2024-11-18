
library(tidyverse)
library(readxl)

m2h <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/mouse.txt", col_names = c("gene", "ortholog"))
h2m <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/human.txt", col_names = c("gene", "ortholog"))

gsets <- c()

s1 <- read_xlsx("41586_2017_BFnature24297_MOESM3_ESM.xlsx", sheet = 1)
s1$foldchange <- s1$foldchange %>% as.numeric()
#up
genes_hs <- s1 %>% filter(foldchange >= 2) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("hangauer-2017_DTP-up", genes_hs, "human"))
gsets <- rbind(gsets, cbind("hangauer-2017_DTP-up", genes_mm, "mouse"))
#down
genes_hs <- s1 %>% filter(foldchange <= 0.5) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("hangauer-2017_DTP-down", genes_hs, "human"))
gsets <- rbind(gsets, cbind("hangauer-2017_DTP-down", genes_mm, "mouse"))

gsets_tib <- tibble(gs_name = gsets[,1],
                    gene = gsets[,2],
                    species = gsets[,3])
write_tsv(x = gsets_tib, file = "markers.txt")
