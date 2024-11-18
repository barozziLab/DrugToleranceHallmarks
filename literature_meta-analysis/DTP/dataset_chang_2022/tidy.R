
library(tidyverse)
library(readxl)

m2h <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/mouse.txt", col_names = c("gene", "ortholog"))
h2m <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/human.txt", col_names = c("gene", "ortholog"))

gsets <- c()

#Table S1: Mesenchymal DTP DEGs
s1 <- read_xlsx("cd-20-1265_supplementary_tables_suppst1.xlsx", sheet = 1, skip = 2)
#up
genes_hs <- s1 %>% filter(`Fold Change` > 0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_DTP_Mesenchymal_all-up", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_DTP_Mesenchymal_all-up", genes_mm, "mouse"))
#down
genes_hs <- s1 %>% filter(`Fold Change` < 0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_DTP_Mesenchymal_all-down", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_DTP_Mesenchymal_all-down", genes_mm, "mouse"))

#Table S2: Luminal-DTP DEGs
s2 <- read_xlsx("cd-20-1265_supplementary_tables_suppst1.xlsx", sheet = 2, skip = 2)
#up
genes_hs <- s2 %>% filter(`Fold Change` > 0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_DTP_Luminal_all-up", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_DTP_Luminal_all-up", genes_mm, "mouse"))
#down
genes_hs <- s2 %>% filter(`Fold Change` < 0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_DTP_Luminal_all-down", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_DTP_Luminal_all-down", genes_mm, "mouse"))

#Table S7: Mesenchymal DTP Unique DEGs
s7 <- read_xlsx("cd-20-1265_supplementary_tables_suppst1.xlsx", sheet = 7, skip = 2)
#up
genes_hs <- s7 %>% filter(Direction == "+") %>% pull(Genes) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_DTP_Mesenchymal_unique-up", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_DTP_Mesenchymal_unique-up", genes_mm, "mouse"))
#down
genes_hs <- s7 %>% filter(Direction == "-") %>% pull(Genes) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_DTP_Mesenchymal_unique-down", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_DTP_Mesenchymal_unique-down", genes_mm, "mouse"))

#Table S8: Luminal DTP Unique DEGs
s8 <- read_xlsx("cd-20-1265_supplementary_tables_suppst1.xlsx", sheet = 8, skip = 2)
#up
genes_hs <- s8 %>% filter(Direction == "+") %>% pull(Genes) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_DTP_Luminal_unique-up", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_DTP_Luminal_unique-up", genes_mm, "mouse"))
#down
genes_hs <- s8 %>% filter(Direction == "-") %>% pull(Genes) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_DTP_Luminal_unique-down", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_DTP_Luminal_unique-down", genes_mm, "mouse"))

#Table S9: DEGs between NPY1R-high and NPY1R-low
s9 <- read_xlsx("cd-20-1265_supplementary_tables_suppst1.xlsx", sheet = 9, skip = 2)
#up
genes_hs <- s9 %>% filter(`Fold change` > 0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_pre-DTP_NPY1R-high-up", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_pre-DTP_NPY1R-high-up", genes_mm, "mouse"))
#down
genes_hs <- s9 %>% filter(`Fold change` < 0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_pre-DTP_NPY1R-high-down", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_pre-DTP_NPY1R-high-down", genes_mm, "mouse"))

#Genes in 5a
genes_5a <- read_tsv("fig_5a.txt", col_names = c("gene"))
genes_hs <- genes_5a %>% pull(gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("chang-2022_pre-DTP-up", genes_hs, "human"))
gsets <- rbind(gsets, cbind("chang-2022_pre-DTP-up", genes_mm, "mouse"))

gsets_tib <- tibble(gs_name = gsets[,1],
                    gene = gsets[,2],
                    species = gsets[,3])
write_tsv(x = gsets_tib, file = "markers.txt")
