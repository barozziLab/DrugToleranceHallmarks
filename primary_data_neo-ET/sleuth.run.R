## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: sleuth.run.R
##
## Description: 
#~
## Sleuth analysis on paired, primary
## bulk RNA-seq data
##
## Authors: 
#~
## Stephan Gruener & Iros Barozzi
##
## License: 
#~
## GNU GPL v3
## Copyright 2024 
## Copyright Iros Barozzi Stephan Gruener
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Notes:
#~
##
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


source("sleuth.libs.R")

source("sleuth.gene_sets.R")

col_pal <- list(treatment = c("pre" = "darkgrey", "post" = "purple"),
                response = c("Responder" = "darkgreen", "Non_responder" = "orange", "Stable" = "lightblue", "NA" = "lightgrey"))


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


#~~~~~~~~~~~~~~~~~~
#Sleuth preparation
###################

outFolder <- "sleuth.results/"
cmd <- paste("mkdir ", outFolder, sep="")
system(cmd, intern=FALSE)

s2c <- data.frame(sample=d_meta_sample$sample_id,
                  patient=d_meta_sample$patient_id,
                  treatment=d_meta_sample$treatment,
                  response=d_meta_sample$response,
                  path=d_meta_sample$path, stringsAsFactors=FALSE)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE, target_mapping = ttg, aggregation_column = "gene_symbol")

# PCA (Sleuth)

outF <- paste(outFolder, "/sleuth.pca-internal.treatment.pdf", sep="")
pdf(outF, width=5.5, height=4)
p <- plot_pca(so, use_filtered = TRUE, units = 'est_counts', color_by = 'treatment') + theme_bw()
plot(p)
dev.off()

outF <- paste(outFolder, "/sleuth.pca-internal.response.pdf", sep="")
pdf(outF, width=6.5, height=4)
p <- plot_pca(so, use_filtered = TRUE, units = 'est_counts', color_by = 'response') + theme_bw()
plot(p)
dev.off()

# Clustering and PCA (using base R) of samples [at least one sample w/TPM >= 1]

exp_m <- sleuth_to_matrix(so, "obs_norm", "tpm")
exp_m <- exp_m[apply(exp_m, 1 ,max) >=1,]
col_names <- d_meta_sample %>% arrange(sample_id) %>% mutate(short_id = paste(patient_id, treatment, sep = "_")) %>% pull(short_id)
colnames(exp_m) <- col_names

exp_m_log <- log2(exp_m + 0.01)
outF <- paste(outFolder, "/sleuth.results.expr_matrix.log2.rds", sep="")
saveRDS(round(exp_m_log, 3) %>% as.data.frame(), file = outF)

#hclust

dst <- as.dist(1 - cor(exp_m, method="spearman"))
hc <- hclust(dst, method="ward.D2")

outF <- paste(outFolder, "/sleuth.hc-r.pdf", sep="")
pdf(outF, width=10, height=5)
plot(hc)
dev.off()

#pca

res.pca <- prcomp(t(exp_m_log), center = TRUE, scale = TRUE)

screen_plot <- fviz_eig(res.pca)
outF <- paste(outFolder, "/sleuth.pca-r.scree.pdf", sep="")
pdf(outF, width = 4, height = 3)
plot(screen_plot)
dev.off()

res_pca_x <- res.pca$x %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "id") %>% 
  as_tibble()

d_meta_sample_tmp <- d_meta_sample %>% 
  rowwise() %>% 
  mutate(id = paste0(patient_id, "_", treatment)) %>% 
  select(-c("sample_id"))

res_pca_x <- res_pca_x %>%
  left_join(d_meta_sample_tmp, by = "id")

top_pcs <- 5
pca_plots <- list()
for (i in 1:top_pcs) {
  if (i != 2) {
    pc_y <- paste0("PC", i)
    pca_plots[[as.character(i)]][["treatment"]] <- res_pca_x %>% 
      ggplot(aes(x=PC2, y=get(pc_y), col=treatment)) + 
      geom_point() + 
      ylab(pc_y) +
      scale_color_manual(values = as.vector(col_pal$treatment)) +
      theme_bw()
    pca_plots[[as.character(i)]][["response"]] <- res_pca_x %>% 
      ggplot(aes(x=PC2, y=get(pc_y), col=response)) + 
      geom_point() + 
      ylab(pc_y) +
      scale_color_manual(values = as.vector(col_pal$response)) +
      theme_bw()
  }
}

outF <- paste(outFolder, "/sleuth.pca-r.scatter.treatment.pdf", sep="")
pdf(outF, width = 5, height = 4)
for (i in 1:top_pcs) {
  if (i != 2) {
    pc_y <- paste0("PC", i)
    plot(pca_plots[[as.character(i)]][["treatment"]])
  }
}
dev.off()

outF <- paste(outFolder, "/sleuth.pca-r.scatter.response.pdf", sep="")
pdf(outF, width = 5.5, height = 4)
for (i in 1:top_pcs) {
  if (i != 2) {
    pc_y <- paste0("PC", i)
    plot(pca_plots[[as.character(i)]][["response"]])
  }
}
dev.off()

res_pca_x_long <- res_pca_x %>% 
  select(id, treatment, response, patient_id, PC1, PC2, PC3, PC4, PC5) %>% 
  pivot_longer(cols = -c(id, treatment, response, patient_id), names_to = "PC")

pca_boxplots <- list()

pca_boxplots$treat <- res_pca_x_long %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot(outlier.size = 0) +
  facet_grid(~PC) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.1) +
  ylim(-100, 100) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pca_boxplots$resp <- res_pca_x_long %>%
  dplyr::filter(response %in% c("Responder", "Non_responder")) %>%
  ggplot(aes(x=response, y=value, fill=response)) +
  geom_boxplot(outlier.size = 0) +
  facet_grid(~PC) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.1) +
  ylim(-100, 100) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

pca_boxplots$both <- res_pca_x_long %>%
  dplyr::filter(response %in% c("Responder", "Non_responder")) %>%
  ggplot(aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot(outlier.size = 0) +
  facet_grid(response~PC) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.1) +
  ylim(-100, 100) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

outF <- paste(outFolder, "/sleuth.pca-r.boxplot.treatment.pdf", sep="")
pdf(outF, width = 5, height = 3)
plot(pca_boxplots$treat)
dev.off()

outF <- paste(outFolder, "/sleuth.pca-r.boxplot.response.pdf", sep="")
pdf(outF, width = 5, height = 3.5)
plot(pca_boxplots$resp)
dev.off()

outF <- paste(outFolder, "/sleuth.pca-r.boxplot.both.pdf", sep="")
pdf(outF, width = 5, height = 4)
plot(pca_boxplots$both)
dev.off()

#Mann-Whitney Tests

ps_treat <- res_pca_x_long %>% 
  group_by(PC) %>% 
  summarize(treatment_mw_p = wilcox.test(value ~ treatment)$p.value)

ps_resp <- res_pca_x_long %>% 
  dplyr::filter(response %in% c("Responder", "Non_responder")) %>%
  group_by(PC) %>% 
  summarize(response_mw_p = wilcox.test(value ~ response)$p.value)

ps <- ps_treat %>%
  left_join(ps_resp, by = "PC")

outF <- paste(outFolder, "/sleuth.pca-r.stats.txt", sep="")
write_tsv(x = ps, file = outF)


#~~~~~~~~~~~~~~~~~~
#Model: Post vs Pre
###################

# To approach the data in a paired way (post vs pre) 
# Use patient as a null model to exclude differences between them

so <- sleuth_fit(so, ~patient, 'patient')
so <- sleuth_fit(so, ~patient + treatment, 'full')

so <- sleuth_lrt(so, 'patient', 'full')
sleuth_table <- sleuth_results(obj = so, test = 'patient:full', test_type = 'lrt', show_all = FALSE)

#sleuth_table_full <- sleuth_results(obj = so, test = 'patient:full', test_type = 'lrt', show_all = FALSE, pval_aggregate = FALSE)
#plot_bootstrap(so, "ENST00000396123.2", units = "est_counts", color_by = "treatment")

# Calculate overall fold-changes without pairing

so <- sleuth_wt(so, 'treatmentpre', 'full')
sleuth_table_wt <- sleuth_results(obj = so, test = 'treatmentpre', test_type = 'wt', show_all = FALSE, pval_aggregate = FALSE)
sleuth_b <- sleuth_table_wt %>% group_by(gene_symbol) %>% summarise(b = -mean(b)) %>% rename(target_id = gene_symbol)

# And add to the table, along with PA annotations

sleuth_table <- sleuth_table %>% as_tibble %>% 
                  left_join(sleuth_b, by = "target_id")

sleuth_table <- sleuth_table %>% 
                  mutate(pa_up = target_id %in% gs_lst$pa$PA_UP) %>% 
                  mutate(pa_down = target_id %in% gs_lst$pa$PA_DOWN)

outF <- paste(outFolder, "/sleuth.results.txt", sep="")
write_tsv(sleuth_table, file = outF)

# Save list of DEGs [arbitrary threshold on b]

degs <- list()
degs[["up"]] <- sleuth_table %>% filter(qval <= 0.05 & b >= 0.5) %>% pull(target_id) %>% sort()
degs[["down"]] <- sleuth_table %>% filter(qval <= 0.05 & b <= -0.5) %>% pull(target_id) %>% sort()
for (cl in names(degs)) {
  outF <- paste(outFolder, "/sleuth.results.degs.", cl, ".txt", sep="")
  write.table(degs[[cl]], file = outF, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Model: Responder vs Non-responder
##################################

#exclude response [NA, Stable]

d_meta_sample_filt <- d_meta_sample %>% filter(!response %in% c(NA, "Stable"))

s2c_filt <- data.frame(sample=d_meta_sample_filt$sample_id,
                       patient=d_meta_sample_filt$patient_id,
                       treatment=d_meta_sample_filt$treatment,
                       response=d_meta_sample_filt$response,
                       path=d_meta_sample_filt$path, stringsAsFactors=FALSE)

so_filt <- sleuth_prep(s2c_filt, extra_bootstrap_summary = TRUE, target_mapping = ttg, aggregation_column = "gene_symbol")

#note: can't do paired [system is computationally singular]
#note: excluding patient factor

so_filt <- sleuth_fit(so_filt, ~ treatment, 'treat')

#compare treat+response to treat model
so_filt <- sleuth_fit(so_filt, ~treatment + response, 'treat_resp')
so_filt <- sleuth_lrt(so_filt, 'treat', 'treat_resp')

sleuth_table_resp <- sleuth_results(obj = so_filt, test = 'treat:treat_resp', test_type = 'lrt', show_all = FALSE)

# Calculate overall fold-changes without pairing

so_filt <- sleuth_wt(obj = so_filt, which_beta = 'responseResponder', which_model = 'treat_resp')
sleuth_table_resp_wt <- sleuth_results(obj = so_filt, test = 'responseResponder', test_type = 'wt', which_model = 'treat_resp', show_all = FALSE, pval_aggregate = FALSE)
sleuth_b <- sleuth_table_resp_wt %>% group_by(gene_symbol) %>% summarise(b = -mean(b)) %>% rename(target_id = gene_symbol)

# And add to the table, along with PA annotations

sleuth_table_resp <- sleuth_table_resp %>% as_tibble %>% 
  left_join(sleuth_b, by = "target_id")

sleuth_table_resp <- sleuth_table_resp %>% 
  mutate(pa_up = target_id %in% gs_lst$pa$PA_UP) %>% 
  mutate(pa_down = target_id %in% gs_lst$pa$PA_DOWN)

outF <- paste(outFolder, "/sleuth.results.nonresp_vs_resp.txt", sep="")
write_tsv(sleuth_table_resp, file = outF)

