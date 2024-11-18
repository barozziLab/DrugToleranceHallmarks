
library(tidyverse)
library(readxl)

data <- read_xlsx("41586_2023_6252_MOESM1_ESM.signatures.bad_formatting.xlsx")
data <- data %>%
  mutate(gene1 = sapply(data$genes, function(x){strsplit(x," ")[[1]][1]})) %>%
  mutate(gene2 = sapply(data$genes, function(x){strsplit(x," ")[[1]][5]}))

d1 <- data %>% 
  select(table, signature1, gene1) %>% 
  dplyr::filter(!is.na(gene1)) %>% 
  dplyr::rename(gene = gene1, signature = signature1)

d2 <- data %>% 
  select(table, signature2, gene2) %>% 
  dplyr::filter(!is.na(gene2)) %>% 
  dplyr::rename(gene = gene2, signature = signature2)

data_clean <- rbind(d1, d2) %>%
  arrange(table, signature)

write_tsv(x = data_clean, file = "41586_2023_6252_MOESM1_ESM.signatures.good_formatting.txt")

