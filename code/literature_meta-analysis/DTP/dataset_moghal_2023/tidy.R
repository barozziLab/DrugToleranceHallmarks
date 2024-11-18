
library(tidyverse)
library(readxl)

m2h <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/mouse.txt", col_names = c("gene", "ortholog"))
h2m <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/human.txt", col_names = c("gene", "ortholog"))

gsets <- c()

#Supplementary Table 4. Differentially expressed genes (DEGs) between DTP and BL tumors by gene expression microarray.
s4 <- read_xlsx("1-s2.0-S1556086422019645-mmc3.xlsx", sheet = 3)
#up
genes_hs <- s4 %>% filter(log2FC >= 1) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("moghal-2023_DTP-up_microarray", genes_hs, "human"))
gsets <- rbind(gsets, cbind("moghal-2023_DTP-up_microarray", genes_mm, "mouse"))
#down
genes_hs <- s4 %>% filter(log2FC <= -1) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("moghal-2023_DTP-down_microarray", genes_hs, "human"))
gsets <- rbind(gsets, cbind("moghal-2023_DTP-down_microarray", genes_mm, "mouse"))

#Supplementary Table 6. Differentially expressed genes (DEGs) between DTP and BL tumor cells by scRNA-seq.
s6 <- read_xlsx("1-s2.0-S1556086422019645-mmc3.xlsx", sheet = 5)
#up
genes_hs <- s6 %>% filter(avg_lnFC >= 0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("moghal-2023_DTP-up_scrnaseq", genes_hs, "human"))
gsets <- rbind(gsets, cbind("moghal-2023_DTP-up_scrnaseq", genes_mm, "mouse"))
#down
genes_hs <- s6 %>% filter(avg_lnFC <= -0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("moghal-2023_DTP-down_scrnaseq", genes_hs, "human"))
gsets <- rbind(gsets, cbind("moghal-2023_DTP-down_scrnaseq", genes_mm, "mouse"))

#Supplementary Table 7. Differentially expressed genes (DEGs) between DTP-like (DTPL) and BL tumor cells by scRNA-seq.
s7 <- read_xlsx("1-s2.0-S1556086422019645-mmc3.xlsx", sheet = 6)
#up
genes_hs <- s7 %>% filter(avg_lnFC >= 0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("moghal-2023_DTPL-up_scrnaseq", genes_hs, "human"))
gsets <- rbind(gsets, cbind("moghal-2023_DTPL-up_scrnaseq", genes_mm, "mouse"))
#down
genes_hs <- s7 %>% filter(avg_lnFC <= -0) %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("moghal-2023_DTPL-down_scrnaseq", genes_hs, "human"))
gsets <- rbind(gsets, cbind("moghal-2023_DTPL-down_scrnaseq", genes_mm, "mouse"))

#Supplementary Table 8. DTP DEGs common to microarray and scRNA-seq.
s8 <- read_xlsx("1-s2.0-S1556086422019645-mmc3.xlsx", sheet = 7)
#up
genes_hs <- s8 %>% filter(Direction == "UP") %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("moghal-2023_DTP-up_common", genes_hs, "human"))
gsets <- rbind(gsets, cbind("moghal-2023_DTP-up_common", genes_mm, "mouse"))
#down
genes_hs <- s8 %>% filter(Direction == "DOWN") %>% pull(Gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("moghal-2023_DTP-down_common", genes_hs, "human"))
gsets <- rbind(gsets, cbind("moghal-2023_DTP-down_common", genes_mm, "mouse"))

gsets_tib <- tibble(gs_name = gsets[,1],
                    gene = gsets[,2],
                    species = gsets[,3])
write_tsv(x = gsets_tib, file = "markers.txt")
