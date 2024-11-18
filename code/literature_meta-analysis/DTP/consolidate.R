
library(tidyverse)
library(readxl)

m2h <- read_tsv("../../homology_map/20200307_ensembl/mouse.txt", col_names = c("gene", "ortholog"))
h2m <- read_tsv("../../homology_map/20200307_ensembl/human.txt", col_names = c("gene", "ortholog"))

gsets <- c()

#~~~~~~~~~~~
#Latest ones
############

sets <- c("hangauer_2017", "moghal_2023", "sammut_2022", "zhao_2020", "chang_2022")

d_set <- list()
for (set in sets) {
  path <- paste0("dataset_", set, "/markers.txt")
  d_set[[set]] <- read_tsv(path)
}

gsets <- rbind(gsets, do.call(rbind, d_set))

############

#~~~~~~~
#PA Hong
########

#Hong et al. 2019

pa_gsets <- c()

#PA - up
inF <- "../../data/PA_UP.txt"
genes_hs <- read_tsv(inF, col_names = FALSE) %>% pull(1) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
pa_gsets <- rbind(pa_gsets, cbind("hong-2019_preDTP-up", genes_hs, "human"))
pa_gsets <- rbind(pa_gsets, cbind("hong-2019_preDTP-up", genes_mm, "mouse"))

#PA - down
inF <- "../../data/PA_DOWN.txt"
genes_hs <- read_tsv(inF, col_names = FALSE) %>% pull(1) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
pa_gsets <- rbind(pa_gsets, cbind("hong-2019_preDTP-down", genes_hs, "human"))
pa_gsets <- rbind(pa_gsets, cbind("hong-2019_preDTP-down", genes_mm, "mouse"))

pa_gsets_tib <- tibble(gs_name = pa_gsets[,1],
                       gene = pa_gsets[,2],
                       species = pa_gsets[,3])

gsets <- rbind(gsets, pa_gsets_tib)

########

#~~~~
#CTCs
#####

ctc_gsets <- c()

## Yu et al. 2013

#BCa up
inF <- "dataset_yu_2013/DEGs.BCa.EPCAM_vs_IgG.up.txt"
genes_hs <- read_tsv(inF, col_names = FALSE) %>% pull(1) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
ctc_gsets <- rbind(ctc_gsets, cbind("yu-2013_CTC-up", genes_hs, "human"))
ctc_gsets <- rbind(ctc_gsets, cbind("yu-2013_CTC-up", genes_mm, "mouse"))

#BCa down
inF <- "dataset_yu_2013/DEGs.BCa.EPCAM_vs_IgG.down.txt"
genes_hs <- read_tsv(inF, col_names = FALSE) %>% pull(1) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
ctc_gsets <- rbind(ctc_gsets, cbind("yu-2013_CTC-down", genes_hs, "human"))
ctc_gsets <- rbind(ctc_gsets, cbind("yu-2013_CTC-down", genes_mm, "mouse"))

## Aceto et al. 2014

#Clusters vs single-cell up
inF <- "dataset_aceto_2014/DEGs.CL_vs_SC.up.txt"
genes_hs <- read_tsv(inF, col_names = TRUE) %>% pull(1) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
ctc_gsets <- rbind(ctc_gsets, cbind("aceto-2014_CTC-clusters-up", genes_hs, "human"))
ctc_gsets <- rbind(ctc_gsets, cbind("aceto-2014_CTC-clusters-up", genes_mm, "mouse"))

#Clusters vs single-cell down
inF <- "dataset_aceto_2014/DEGs.SC_vs_CL.up.txt"
genes_hs <- read_tsv(inF, col_names = TRUE) %>% pull(1) %>% sort()
genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
ctc_gsets <- rbind(ctc_gsets, cbind("aceto-2014_CTC-clusters-down", genes_hs, "human"))
ctc_gsets <- rbind(ctc_gsets, cbind("aceto-2014_CTC-clusters-down", genes_mm, "mouse"))

ctc_gsets_tib <- tibble(gs_name = ctc_gsets[,1],
                        gene = ctc_gsets[,2],
                        species = ctc_gsets[,3])

gsets <- rbind(gsets, ctc_gsets_tib)

#####

#~~~~~~~~
#Diapause
#########

diapause_gsets <- c()

#Dhimolea et al. 2021

d_dhimolea <- readRDS("dataset_dhimolea_2021/deseq_results.rds")

names_oi <-  c("HCI002_Org-docetaxel_VS_HCI002_Org-vehicle",
               "HCI003_Org-docetaxel_VS_HCI003_Org-vehicle",
               "HCI003_PDX-afatinib_VS_HCI003_PDX-vehicle",     
               "MDAMB231_Org-docetaxel_VS_MDAMB231_Org-vehicle",
               "PM154_Org-docetaxel_VS_PM154_Org-vehicle",
               "PrCa1_Org-docetaxel_VS_PrCa1_Org_A-vehicle",
               "PrCa1_Org-docetaxel_VS_PrCa1_Org_B-vehicle",
               "PrCa1_Xeno-vinblastine_VS_PrCa1_Xeno-vehicle",
               "PrCa3_Org-docetaxel_VS_PrCa3_Org-vehicle",
               "PrCa7_Org-docetaxel_VS_PrCa7_Org-vehicle")

names_oi_short <-  c("Org-HCI002-Docetaxel",
                     "Org-HCI003-Docetaxel",
                     "PDX-HCI003-Afatinib",
                     "Org-MDAMB231-Docetaxel",
                     "Org-PM154-Docetaxel",
                     "Org-PrCa1-Docetaxel",
                     "Org-PrCa1-Docetaxel",
                     "Xeno-PrCa1-Vinblastine",
                     "Org-PrCa3-Docetaxel",
                     "Org-PrCa7-Docetaxel")

for (i in 1:length(names_oi)) {
  name_oi <- names_oi[i]
  #up
  genes_hs <- d_dhimolea[[name_oi]] %>% filter(log2FoldChange >= 1 & padj <= 0.05) %>% arrange(Symbol) %>% pull(Symbol)
  if (length(genes_hs) > 0) {
    genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
    diapause_gsets <- rbind(diapause_gsets, cbind(paste0("dhimolea-2021_Diapause-", names_oi_short[i], "-up"), genes_hs, "human"))
    diapause_gsets <- rbind(diapause_gsets, cbind(paste0("dhimolea-2021_Diapause-", names_oi_short[i], "-up"), genes_mm, "mouse"))
  }
  #down
  genes_hs <- d_dhimolea[[name_oi]] %>% filter(log2FoldChange <= -1 & padj <= 0.05) %>% arrange(Symbol) %>% pull(Symbol)
  if (length(genes_hs) > 0) {
    genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
    diapause_gsets <- rbind(diapause_gsets, cbind(paste0("dhimolea-2021_Diapause-", names_oi_short[i], "-down"), genes_hs, "human"))
    diapause_gsets <- rbind(diapause_gsets, cbind(paste0("dhimolea-2021_Diapause-", names_oi_short[i], "-down"), genes_mm, "mouse"))
  }
}

## Rehman et al. 2021

inFolder <- "dataset_rehman_2021/rehman_et_al_2021/"

inFiles <- c("DTP.down.txt", "DTP.up.txt")
suffixes <- c("down", "up")

for (i in 1:length(inFiles)) {
  inF <- paste(inFolder, inFiles[i], sep = "")
  genes_hs <- read_tsv(inF, col_names = FALSE) %>% pull(1) %>% sort()
  genes_mm <- h2m %>% filter(gene %in% genes_hs) %>% pull(ortholog)
  diapause_gsets <- rbind(diapause_gsets, cbind(paste0("rehman-2021_Diapause-", suffixes[i]), genes_hs, "human"))
  diapause_gsets <- rbind(diapause_gsets, cbind(paste0("rehman-2021_Diapause-", suffixes[i]), genes_mm, "mouse"))
}

diapause_gsets_tib <- tibble(gs_name = diapause_gsets[,1],
                             gene = diapause_gsets[,2],
                             species = diapause_gsets[,3])

gsets <- rbind(gsets, diapause_gsets_tib)

#########

write_tsv(x = gsets, file = "DTPs.markers.txt")

