
library("tidyverse")
library("DESeq2")

#~~~~~~~~~~~~~~~~
#Data Preparation
#################

exp_data <- read_tsv("GSE197780_DARANA_GE_table.txt.gz", skip = 1)

#Gene symbol -> mean across ensembl ids

exp_data_long <- exp_data %>% 
  pivot_longer(cols = -c("ensembl_gene_id", "gene_id")) %>%
  group_by(gene_id, name) %>%
  summarise(value = mean(value)) %>%
  ungroup()

# to wide

exp_data_wide <- exp_data_long %>% 
  pivot_wider(names_from = name, values_from = value)

exp_data_mat <- exp_data_wide %>%
  column_to_rownames(var = "gene_id")

# to integers

exp_data_mat <- round(2^exp_data_mat-1)

# column-order matched metadata

meta <- tibble(
  individual = as.vector(sapply(colnames(exp_data_mat), function(x){strsplit(x, "_")[[1]][1]})),
  treatment = as.vector(sapply(colnames(exp_data_mat), function(x){strsplit(x, "_")[[1]][2]}))
)

meta <- meta %>% 
  mutate(individual = factor(individual)) %>% 
  mutate(treatment = factor(treatment, levels = c("pre", "post")))


#~~~~~~
#DESeq2
#######

se <- SummarizedExperiment(assays = list(counts = as.matrix(exp_data_mat)), 
                           colData = as.data.frame(meta))

dds <- DESeqDataSet(se, design = ~ individual + treatment)

#normalisation
dds <- estimateSizeFactors(dds)

#colData(dds)

#save dds
#out_file <- paste(out_folder, "/dds.rds", sep = "")
#saveRDS(object=dds, file=out_file)

#VST
vsd <- vst(dds, blind = TRUE)

#PCA
pca_out <- plotPCA(vsd, intgroup = c("treatment"))

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

#add treatment to exp_data_long
exp_data_long$individual <- as.vector(sapply(exp_data_long$name, function(x){strsplit(x, "_")[[1]][1]}))
exp_data_long$treatment <- as.vector(sapply(exp_data_long$name, function(x){strsplit(x, "_")[[1]][2]}))

#exp_data_long %>% dplyr::filter(gene_id == "JUN") %>% ggplot(aes(x = treatment, y = value)) + geom_boxplot()
#exp_data_long %>% dplyr::filter(gene_id == "KLK3") %>% ggplot(aes(x = treatment, y = value)) + geom_boxplot()

