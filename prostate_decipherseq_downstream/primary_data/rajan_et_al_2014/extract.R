
library("tidyverse")
library("readxl")

data <- read_xlsx("12885_2014_5118_MOESM3_ESM.xlsx")

data <- data %>% 
  mutate(logFC = as.numeric(logFC)) %>%
  mutate(PValue = as.numeric(PValue)) %>%
  mutate(FDR = as.numeric(FDR))

up_reg <- data %>% 
  dplyr::filter(FDR <= 0.05 & logFC >= log2(1.5)) %>% 
  pull(Symbol) %>% 
  sort()

down_reg <- data %>% 
  dplyr::filter(FDR <= 0.05 & logFC <= -log2(1.5)) %>% 
  pull(Symbol) %>% 
  sort()

write.table(x = up_reg, file = "tax_pre_vs_post.up.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(x = down_reg, file = "tax_pre_vs_post.down.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)

