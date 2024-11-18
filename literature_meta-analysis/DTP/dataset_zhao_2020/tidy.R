
library(tidyverse)
library(readxl)

m2h <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/mouse.txt", col_names = c("gene", "ortholog"))
h2m <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/human.txt", col_names = c("gene", "ortholog"))

gsets <- c()

s1_GSE153183_up <- read_xlsx("Table 1.XLSX", sheet = 1)
genes_hs <- s1_GSE153183_up %>% pull(gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("zhao-2020_DTP-up_GSE153183", genes_hs, "human"))
gsets <- rbind(gsets, cbind("zhao-2020_DTP-up_GSE153183", genes_mm, "mouse"))

s1_GSE153183_down <- read_xlsx("Table 1.XLSX", sheet = 2)
genes_hs <- s1_GSE153183_down %>% pull(gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("zhao-2020_DTP-down_GSE153183", genes_hs, "human"))
gsets <- rbind(gsets, cbind("zhao-2020_DTP-down_GSE153183", genes_mm, "mouse"))

s2_GSE114647_up <- read_xlsx("Table 2.XLSX", sheet = 1)
genes_hs <- s2_GSE114647_up %>% pull(symbol) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("zhao-2020_DTP-up_GSE114647", genes_hs, "human"))
gsets <- rbind(gsets, cbind("zhao-2020_DTP-up_GSE114647", genes_mm, "mouse"))

s2_GSE114647_down <- read_xlsx("Table 2.XLSX", sheet = 2)
genes_hs <- s2_GSE114647_down %>% pull(symbol) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("zhao-2020_DTP-down_GSE114647", genes_hs, "human"))
gsets <- rbind(gsets, cbind("zhao-2020_DTP-down_GSE114647", genes_mm, "mouse"))

#s1/s2 intersection
genes_hs <- sort(intersect(s1_GSE153183_up %>% pull(gene), s2_GSE114647_up %>% pull(symbol)))
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("zhao-2020_DTP-up_both", genes_hs, "human"))
gsets <- rbind(gsets, cbind("zhao-2020_DTP-up_both", genes_mm, "mouse"))
genes_hs <- sort(intersect(s1_GSE153183_down %>% pull(gene), s2_GSE114647_down %>% pull(symbol)))
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("zhao-2020_DTP-down_both", genes_hs, "human"))
gsets <- rbind(gsets, cbind("zhao-2020_DTP-down_both", genes_mm, "mouse"))

s3 <- read_xlsx("Table 3.XLSX", sheet = 1)
genes_hs <- s3 %>% pull(gene) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
gsets <- rbind(gsets, cbind("zhao-2020_DTP-univariateCOX", genes_hs, "human"))
gsets <- rbind(gsets, cbind("zhao-2020_DTP-univariateCOX", genes_mm, "mouse"))

gsets_tib <- tibble(gs_name = gsets[,1],
                    gene = gsets[,2],
                    species = gsets[,3])
write_tsv(x = gsets_tib, file = "markers.txt")
