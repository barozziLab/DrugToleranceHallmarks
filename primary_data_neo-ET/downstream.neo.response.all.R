
## What: this tests for differences between treatment-naive responder vs non-responder (stable + non-responder)

source("downstream.libs.R")

outFolder <- "downstream.neo.response.results/"
cmd <- paste("mkdir", outFolder)
system(cmd)


#~~~~~~~~~~~~~~~~~~~~
#Metadata preparation
#####################

# Create a table to map target_id to gene symbol, so this can be used to aggregate the analysis to gene symbols

d_meta_genes <- read_tsv("TranscriptToGene.txt")
target_id <- d_meta_genes %>% pull(4)
gene_symbol <- d_meta_genes %>% pull(3)
ttg <- data.frame(target_id, gene_symbol = gene_symbol, stringsAsFactors = FALSE)

# Clinical and Samples Metadata Preparation

d_meta_patient <- read_tsv("metadata_clinical.txt")
d_meta_sample <- read_tsv("metadata_samples.txt")
d_meta_sample <- d_meta_sample %>% mutate(path = paste("pseudoAlignment/", sample_id, sep = ""))

d_meta_temp <- d_meta_patient %>% 
  mutate(response = `Responder/Non_responder/Stable`) %>% 
  dplyr::select(patient_id, response)
d_meta_sample <- d_meta_sample %>% left_join(d_meta_temp, by = "patient_id")


#~~~~~~
#Sleuth
#######

#only pre-; responders vs nr/stable (exclude NA)
#subset metadata

d_meta_sample_pre <- d_meta_sample %>% 
  dplyr::filter(treatment == "pre" & !is.na(response)) %>% 
  mutate(respose_2cls = ifelse(response == "Responder", "Responder", "Non_responder"))

s2c <- data.frame(sample=d_meta_sample_pre$sample_id,
                  patient=d_meta_sample_pre$patient_id,
                  response=d_meta_sample_pre$respose_2cls,
                  path=d_meta_sample_pre$path, stringsAsFactors=FALSE)

outF <- paste(outFolder, "/sleuth.rds", sep="")
if (file.exists(outF)) {
  so <- readRDS(file = outF)
} else {
  so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = ttg, aggregation_column = "gene_symbol")
  saveRDS(object = so, file = outF)
}

# PCA (Sleuth)

p2 <- plot_pca(so, use_filtered = TRUE, units = 'est_counts', color_by = 'response', pc_x = 2, pc_y = 1) + theme_bw()
p3 <- plot_pca(so, use_filtered = TRUE, units = 'est_counts', color_by = 'response', pc_x = 2, pc_y = 3) + theme_bw()
p4 <- plot_pca(so, use_filtered = TRUE, units = 'est_counts', color_by = 'response', pc_x = 2, pc_y = 4) + theme_bw()

outF <- paste(outFolder, "/sleuth.pca-internal.pdf", sep="")
pdf(outF, width=5, height=7)
plot(p2 / p3 / p4)
dev.off()

# Clustering and PCA (using base R) of samples [at least one sample w/TPM >= 1]

exp_m <- sleuth_to_matrix(so, "obs_norm", "tpm")
exp_m <- exp_m[apply(exp_m, 1 ,max) >=1,]
col_names <- d_meta_sample_pre %>% arrange(sample_id) %>% mutate(short_id = paste(patient_id, treatment, respose_2cls, sep = "_")) %>% pull(short_id)
colnames(exp_m) <- col_names

exp_m_log <- log2(exp_m + 0.01)
outF <- paste(outFolder, "/sleuth.results.expr_matrix.log2.rds", sep="")
saveRDS(round(exp_m_log, 3) %>% as.data.frame(), file = outF)

dst <- as.dist(1 - cor(exp_m, method="spearman"))
hc <- hclust(dst, method="ward.D2")

outF <- paste(outFolder, "/sleuth.hc-r.pdf", sep="")
pdf(outF, width=10, height=7)
plot(hc)
dev.off()


#~~~~~~~~~~~~~~
#Models Fitting
###############

so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_fit(so, ~response, 'response')

so <- sleuth_lrt(so, 'reduced', 'response')
sleuth_table <- sleuth_results(obj = so, test = 'reduced:response', test_type = 'lrt', show_all = FALSE)

# Generate expression table collapsed at gene symbol

gene_map <- read_tsv("TranscriptToGene.txt")
colnames(gene_map) <- c("ensembl_gene", "ensembl_tx", "symbol", "ensembl_tx_v")

dset_ens <- exp_m_log %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(ensembl_tx_v = rowname)

dset_ens_long <- dset_ens %>% 
  pivot_longer(!ensembl_tx_v, names_to = "sample", values_to = "expr") %>%
  right_join(gene_map, by = "ensembl_tx_v") %>%
  dplyr::filter(!is.na(sample))

dset_ens_long_avg <- dset_ens_long %>% 
  group_by(sample, symbol) %>% 
  summarise(expr_avg = mean(expr))

dset_tib <- dset_ens_long_avg %>% pivot_wider(names_from = sample, values_from = expr_avg)
dset <- dset_tib %>% 
  select(-c("symbol")) %>% 
  as.data.frame()
rownames(dset) <- dset_tib %>% pull(symbol)

# Calculate fold-changes

dset_long <- dset %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c("gene")) %>%
  rowwise() %>%
  mutate(patient = strsplit(name, split = "_pre_")[[1]][1]) %>%
  mutate(response = strsplit(name, split = "_pre_")[[1]][2])

dset_long_summary <- dset_long %>% 
  group_by(gene, response) %>% 
  summarise(mean = mean(value), n = n()) %>%
  ungroup()

dset_fcs <- dset_long_summary %>%
  select(-c("n")) %>%
  pivot_wider(names_from = "response", values_from = "mean") %>%
  select(gene, Responder, Non_responder) %>%
  mutate(log2fc = Responder - Non_responder)

# And add to the table, then save

sleuth_table <- sleuth_table %>% 
  as_tibble() %>%
  mutate(DEG = qval <= 0.05) %>%
  dplyr::rename(gene = target_id) %>%
  left_join(dset_fcs %>% select(gene, log2fc), by = "gene")

outF <- paste(outFolder, "/sleuth.results.txt", sep="")
write_tsv(sleuth_table, file = outF)

#plotting single genes
#plot_bootstrap(so, "ENST00000498016.1", units = "est_counts", color_by = "response")
#plot_bootstrap(so, "ENST00000631379.1", units = "est_counts", color_by = "response")

###############


#~~~~~~~~~~~~~
#DecoupleR TFs
##############

net_tfs <- get_collectri(organism='human', split_complexes=FALSE)

#sample activities
sample_acts_tfs <- run_ulm(mat=dset, 
                           net=net_tfs, 
                           .source='source', 
                           .target='target', 
                           .mor='mor', 
                           minsize = 10)

sample_acts_mat <- sample_acts_tfs %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

#activities directly from DEGs
#use fold changes from differential expression

degs_fc_mat <- dset_fcs %>% 
  select(gene, log2fc) %>% 
  column_to_rownames("gene")

contrast_acts <- run_ulm(mat=degs_fc_mat, 
                         net=net_tfs, 
                         .source='source', 
                         .target='target', 
                         .mor='mor', 
                         minsize = 10)

#correct for multiple hypotheses testing
contrast_acts <- contrast_acts %>%
  mutate(p_value_BH = p.adjust(contrast_acts$p_value, method = "BH"))

#only significant ones
tfs <- contrast_acts %>% 
  dplyr::filter(p_value_BH <= 0.05) %>%
  pull(source)

#extract the matrix of this TFs for
sample_acts_mat_diff <- sample_acts_mat[,tfs]

#transpose and scale
sample_acts_mat_diff <- t(sample_acts_mat_diff)
cnames <- colnames(sample_acts_mat_diff)
sample_acts_mat_diff <- t(apply(sample_acts_mat_diff, 1, scale))
colnames(sample_acts_mat_diff) <- cnames

response_variable <- sapply(colnames(sample_acts_mat_diff), function(x){strsplit(x, split = "_pre_")[[1]][2]}) %>% as.vector()

ha_cols <- list(response = c("Responder" = "darkgreen", "Non_responder" = "orange"))

ha <- HeatmapAnnotation(response = response_variable,
                        col = ha_cols)

HM <- Heatmap(sample_acts_mat_diff,
              heatmap_legend_param = list(title = "TF Activity (z-score)"),
              col = colorRamp2(c(-2,0,2), c("navy", "white", "red")),
              row_gap = unit(c(3), "mm"),
              row_names_gp = grid::gpar(fontsize = 8),
              column_names_gp = grid::gpar(fontsize = 8),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(3, "cm"),
              border = TRUE,
              top_annotation = ha,
              clustering_distance_rows = "pearson",
              clustering_distance_columns = "euclidean")

outF <- paste(outFolder,"/sleuth.decoupleR.tfs.heatmap.pdf", sep = "")
pdf(outF, width=8, height=9)
draw(HM)
dev.off()

#same as box plots
sample_acts_filt <- sample_acts_mat[,tfs] %>% 
  as.data.frame %>% 
  rownames_to_column(var = "id") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c("id"), names_to = "pathway") %>% 
  rowwise() %>%
  mutate(patient = strsplit(id, split = "_pre_")[[1]][1]) %>%
  mutate(response = strsplit(id, split = "_pre_")[[1]][2])

bplot <- ggplot(data = sample_acts_filt, aes(x=response, y=value)) +
  facet_wrap(~pathway, ncol = 1) +
  geom_boxplot(outlier.size = 0) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.1) +
  coord_flip() +
  theme_bw()

outF <- paste(outFolder,"/decoupleR.tfs.boxplot.pdf", sep = "")
pdf(outF, width=3, height=32)
plot(bplot)
dev.off()

#save a table of full results
outF <- paste(outFolder,"/sleuth.decoupleR.tfs.txt", sep = "")
write_tsv(x = contrast_acts, file = outF)


#~~~~~~~~~~~~~~~~~~~
#DecoupleR Signaling
####################

net_sig <- get_progeny(organism = 'human', top = 500)

#sample activities
sample_acts_sig <- run_ulm(mat=dset, 
                           net=net_sig, 
                           .source='source', 
                           .target='target', 
                           .mor='weight', 
                           minsize = 5)

sample_acts_mat <- sample_acts_sig %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

#activities directly from DEGs
#use b value from differential expression
contrast_acts <- run_ulm(mat=degs_fc_mat, 
                         net=net_sig, 
                         .source='source', 
                         .target='target', 
                         .mor='weight', 
                         minsize = 5)
#correct for multiple hypotheses testing
contrast_acts <- contrast_acts %>%
  mutate(p_value_BH = p.adjust(contrast_acts$p_value, method = "BH"))

#only significant ones, transpose and scale
w <- contrast_acts %>% dplyr::filter(p_value_BH <= 0.05) %>% pull(source)
sample_acts_mat_scaled <- t(sample_acts_mat[,w])
sample_acts_mat_scaled <- t(apply(sample_acts_mat_scaled, 1, scale))
colnames(sample_acts_mat_scaled) <- rownames(sample_acts_mat)

HM <- Heatmap(sample_acts_mat_scaled,
              heatmap_legend_param = list(title = "Sig. Pathway\n Activity (z-score)"),
              col = colorRamp2(c(-2,0,2), c("navy", "white", "red")),
              row_gap = unit(c(3), "mm"),
              row_names_gp = grid::gpar(fontsize = 8),
              column_names_gp = grid::gpar(fontsize = 8),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(3, "cm"),
              border = TRUE,
              top_annotation = ha)

outF <- paste(outFolder,"/decoupleR.signaling.heatmap.pdf", sep = "")
pdf(outF, width=9, height=4.25)
draw(HM)
dev.off()

#same as box plots
sample_acts_filt <- sample_acts_mat[,w] %>% 
  as.data.frame %>% 
  rownames_to_column(var = "id") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c("id"), names_to = "pathway") %>% 
  rowwise() %>%
  mutate(patient = strsplit(id, split = "_pre_")[[1]][1]) %>%
  mutate(response = strsplit(id, split = "_pre_")[[1]][2])

bplot <- ggplot(data = sample_acts_filt, aes(x=response, y=value)) +
  facet_wrap(~pathway, ncol = 1) +
  geom_boxplot(outlier.size = 0) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.1) +
  coord_flip() +
  theme_bw()

outF <- paste(outFolder,"/decoupleR.signaling.boxplot.pdf", sep = "")
pdf(outF, width=3, height=4)
plot(bplot)
dev.off()

#save a table of full results
outF <- paste(outFolder,"/decoupleR.signaling.txt", sep = "")
write_tsv(x = contrast_acts, file = outF)

