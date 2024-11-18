
library(tidyverse)

m2h <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/mouse.txt", col_names = c("gene", "ortholog"))
h2m <- read_tsv("~/Documents/LBNL/datasets/homology_maps/20200307_ensembl/human.txt", col_names = c("gene", "ortholog"))

gsets <- c()

d <- read_tsv("genes.txt")
sets <- d %>% pull(Set) %>% unique()
for (set in sets) {
  genes_hs <- d %>% filter(Set == set) %>% pull(Gene) %>% sort()
  genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
  gset_name <- paste0("sammut-2022_", set)
  gsets <- rbind(gsets, cbind(gset_name, genes_hs, "human"))
  gsets <- rbind(gsets, cbind(gset_name, genes_mm, "mouse"))
}

gsets_tib <- tibble(gs_name = gsets[,1],
                    gene = gsets[,2],
                    species = gsets[,3])
write_tsv(x = gsets_tib, file = "markers.txt")
