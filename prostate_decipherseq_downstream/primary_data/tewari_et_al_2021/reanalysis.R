
library("tidyverse")
library("DESeq2")

#~~~~~~~~~~~~~~~~
#Data Preparation
#################

meta <- read_tsv(file = "1-s2.0-S2211124721011098-mmc2.sheet2.txt")

exp_data <- read_tsv("celrep_109665_mmc9.txt")

# split symbol/ensembl
# to long format
# add info about the individual

strsplit_custom <- function(gene_id) {
  res <- strsplit(gene_id, "_")[[1]]
  return(res[length(res)])
}

exp_data_long <- exp_data %>% 
  rowwise() %>% 
  mutate(gene_symbol = strsplit_custom(gene_id)) %>% 
  pivot_longer(cols = -c("gene_id", "gene_symbol")) %>% 
  rowwise() %>% 
  mutate(sample = gsub("ER_|NR_", "", name)) %>% 
  mutate(sample = strsplit(sample, "_")[[1]][1]) %>% 
  ungroup()

# group by gene and individual
# geometric mean --> exp(mean(log(x+1)))

exp_data_long <- exp_data_long %>% 
  group_by(gene_symbol, sample) %>% 
  summarise(n = n(), value = exp(mean(log(value+1)))-1) %>%
  mutate(value = as.integer(value)) %>%
  ungroup()

# to wide

exp_data_wide <- exp_data_long %>% 
  select(-c("n")) %>% 
  pivot_wider(names_from = sample, values_from = value)

exp_data_mat <- exp_data_wide %>%
  column_to_rownames(var = "gene_symbol")

# column-order matched metadata
# order is already correct, but some samples are missing in the RNA (14-283-64)

meta <- meta %>% 
  dplyr::filter(SID2 %in% colnames(exp_data_wide)) %>%
  mutate(response = factor(response, levels = c("Non", "Exceptional")))


#~~~~~~
#DESeq2
#######

se <- SummarizedExperiment(assays = list(counts = as.matrix(exp_data_mat)), 
                           colData = as.data.frame(meta))

dds <- DESeqDataSet(se, design = ~ response)

#normalisation
dds <- estimateSizeFactors(dds)

#colData(dds)

#save dds
#out_file <- paste(out_folder, "/dds.rds", sep = "")
#saveRDS(object=dds, file=out_file)


#VST
vsd <- vst(dds, blind = TRUE)

#PCA
pca_out <- plotPCA(vsd, intgroup = c("response"))

#alternative plot
#p <- ggplot(pca_out$data %>% mutate(response = as.factor(response)), aes(x = PC1, y = PC2, col = response)) + 
#  geom_point(size = 2) +
#  theme_bw()

pdf("deseq2.results.pca.pdf", width = 4, height = 3.5)

plot(pca_out +
geom_point(size = 2) +
  theme_bw())

dev.off()

#DESeq2
dds <- DESeq(dds)

#results
res_tib <- results(dds) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble()

#add boolean for significance, based on defined thresholds

thresh_lfc <- log2(1.5)
thresh_fdr <- 0.05

res_tib <- res_tib %>% 
  mutate(significant = abs(log2FoldChange) >= thresh_lfc & padj <= thresh_fdr) %>%
  mutate(significant = ifelse(is.na(significant), FALSE, significant)) %>% 
  mutate(significant_up = log2FoldChange >= thresh_lfc & padj <= thresh_fdr) %>%
  mutate(significant_up = ifelse(is.na(significant_up), FALSE, significant_up)) %>% 
  mutate(significant_down = log2FoldChange <= -thresh_lfc & padj <= thresh_fdr) %>%
  mutate(significant_down = ifelse(is.na(significant_down), FALSE, significant_down))

#save
write_tsv(res_tib, file = "deseq2.results.txt")

#add response to exp_data_long
exp_data_long <- exp_data_long %>% 
  left_join(meta %>% select(SID2, response) %>% dplyr::rename(sample = SID2), by = "sample")

#exp_data_long %>% dplyr::filter(gene_symbol == "ANPEP") %>% ggplot(aes(x = response, y = log2(value+1))) + geom_boxplot()
#exp_data_long %>% dplyr::filter(gene_symbol == "TDO2") %>% ggplot(aes(x = response, y = log2(value+1))) + geom_boxplot()

