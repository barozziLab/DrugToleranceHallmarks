
library(tidyverse)

outFolder <- "downstream.pervasiveness.results/"


#~~~~~~~~~~~~~~~~~~~~~
#Gene Sets Preparation
######################

## PA

gs_pa <- list()

in_f <- "../data/PA_UP.txt"
gs_pa[["PA_UP"]] <- read_tsv(in_f, col_names = FALSE) %>% pull(1) %>% sort()

in_f <- "../data/PA_DOWN.txt"
gs_pa[["PA_DOWN"]] <- read_tsv(in_f, col_names = FALSE) %>% pull(1) %>% sort()

gs_pa_df <- rbind(cbind("PA_UP", gs_pa[["PA_UP"]]), cbind("PA_DOWN", gs_pa[["PA_DOWN"]]))
gs_pa_tib <- tibble(gs_name = gs_pa_df[,1], gene_symbol = gs_pa_df[,2])

t2g <- list()
t2g[["Pre-Adapted"]] <- gs_pa_tib

## Breast, DElegate DEGs

degs_res <- readRDS("../breast_seurat/MCF7_cca_kit_5k_feat/deg_results_deseq.rds")

gs_degs <- c()
gs_degs_df <- c()

for (set_name in names(degs_res)) {
  
  set_name_up <- paste0(set_name, "_top1k-up")
  gs_degs[[set_name_up]] <- degs_res[[set_name]] %>% dplyr::filter(padj <= 0.05 & log_fc >= log2(1.5)) %>% head(1000) %>% pull(feature) %>% sort() %>% unique()
  set_name_down <- paste0(set_name, "_top1k-down")
  gs_degs[[set_name_down]] <- degs_res[[set_name]] %>% dplyr::filter(padj <= 0.05 & log_fc <= -log2(1.5)) %>% head(1000) %>% pull(feature) %>% sort() %>% unique()
  
  gs_degs_df <- rbind(gs_degs_df, cbind(set_name_up, gs_degs[[set_name_up]]))
  gs_degs_df <- rbind(gs_degs_df, cbind(set_name_down, gs_degs[[set_name_down]]))
  
}

gs_degs_tib <- tibble(gs_name = gs_degs_df[,1], gene_symbol = gs_degs_df[,2])


#~~~~~~~~~~~~~~~~~~~~~
#STRING-DB Preparation
######################

string_db <- read_delim("../data/20240411_hs.9606.protein.links.v12.0.txt.gz", delim = " ")

#map ENSP to gene_symbol
string_db_ann <- read_tsv("../data/20240411_hs.9606.protein.info.v12.0.txt.gz")
colnames(string_db_ann)[1] <- "protein_id"
#
string_db_gene <- string_db %>% 
  left_join(string_db_ann %>% dplyr::rename(protein1 = protein_id, gene1 = preferred_name) %>% select(-c("protein_size", "annotation")), by = "protein1") %>%
  left_join(string_db_ann %>% dplyr::rename(protein2 = protein_id, gene2 = preferred_name) %>% select(-c("protein_size", "annotation")), by = "protein2") %>%
  select(gene1, gene2, combined_score) %>% 
  unique()


#~~~~~~~~~~~~~~~~~~~~~~~~~
#Pervasiveness Preparation
##########################

perv_sig <- list()

in_F <- paste(outFolder, "/DTP-signatures.txt", sep = "")
perv_sig[["DTP"]] <- read_tsv(in_F)
perv_sig[["DTP"]] <- perv_sig[["DTP"]] %>% dplyr::mutate(n_sig = ifelse(class == "down", -n, n))

in_F <- paste(outFolder, "/NEO-signatures.txt", sep = "")
perv_sig[["NEO"]] <- read_tsv(in_F)
perv_sig[["NEO"]] <- perv_sig[["NEO"]] %>% dplyr::mutate(n_sig = ifelse(class == "down", -n, n))


#~~~~~~~~
#Analysis
#########

d_input <- list()

combined_score_min <- 500

#prepare data (join and filter)

for (name_sig in names(perv_sig)) {
  
  d_input[[name_sig]] <- string_db_gene %>% 
    inner_join(perv_sig[[name_sig]] %>% dplyr::rename(gene1 = gene_symbol), by = "gene1") %>%
    select(gene1, gene2, combined_score, n, n_sig) %>%
    dplyr::rename(gene1_n = n, gene1_n_sig = n_sig) %>%
    inner_join(perv_sig[[name_sig]] %>% dplyr::rename(gene2 = gene_symbol), by = "gene2") %>%
    select(gene1, gene2, combined_score,gene1_n, gene1_n_sig, n, n_sig) %>%
    dplyr::rename(gene2_n = n, gene2_n_sig = n_sig) %>%
    dplyr::filter(combined_score >= combined_score_min)
  
}

#save input tables
for (name_sig in names(perv_sig)) {
  out_path <- paste0(outFolder, "/string.", name_sig, "-signature.annotations.txt")
  write_tsv(x = d_input[[name_sig]], file = out_path)
}

## Is the density of connection in STRING positively correlated w/pervasiveness?

#observed = count /2 (since links are present in both directions)
#total possible = N(N-1)/2

stats <- c()
genes <- c()

for (name_sig in names(perv_sig)) {
  
  for (n_sig in 1:5) {
  
    d_filt <- d_input[[name_sig]] %>% dplyr::filter(abs(gene1_n) >= n_sig & abs(gene2_n) >= n_sig)
    
    obs <- nrow(d_filt) / 2 
    nodes <- c(d_filt %>% pull(gene1), d_filt %>% pull(gene2)) %>% sort() %>% unique()
    N <- length(nodes)
    tot <- N * (N-1) / 2
    
    stats <- rbind(stats, c(name_sig, n_sig, N, obs, tot))
    
    genes <- rbind(genes, cbind(name_sig, n_sig, nodes))
  
  }
  
}

stats <- tibble(name = stats[,1],
                pervasiveness_min = as.numeric(stats[,2]),
                genes_N = as.numeric(stats[,3]),
                links_obs = as.numeric(stats[,4]),
                links_max = as.numeric(stats[,5]))
  
stats <- stats %>% 
  mutate(links_obs_f = links_obs / links_max)

genes <- tibble(name = genes[,1],
                pervasiveness_min = as.numeric(genes[,2]),
                gene = genes[,3])

#save stats
out_path <- paste0(outFolder, "/string.stats.txt")
write_tsv(x = stats, file = out_path)

#save genes
out_path <- paste0(outFolder, "/string.stats.genes.txt")
write_tsv(x = genes, file = out_path)

#plot density ~ pervasiveness
plt <- stats %>% 
  ggplot(aes(x = pervasiveness_min, y = links_obs_f)) +
  geom_line() + 
  facet_wrap(~name) + 
  theme_bw() +
  xlab("Pervasiveness (min)") + 
  ylab("Density of links (STRING-DB)") +
  ylim(0, 0.15)

out_F <- paste(outFolder, "/string.stats.pdf", sep = "")
pdf(out_F, width = 3.5, height = 2.5)
plot(plt)
dev.off()

#Annotations

genes_anno <- genes %>% 
  group_by(gene, name) %>% 
  summarise(pervasiveness_core = max(pervasiveness_min)) %>% 
  pivot_wider(names_from = name, values_from = pervasiveness_core) %>% 
  ungroup() %>%
  mutate(DTP = ifelse(is.na(DTP), 0, DTP)) %>%
  mutate(NEO = ifelse(is.na(NEO), 0, NEO))

#genes_anno %>% dplyr::filter(NEO != 0) %>% ggplot(aes(x = as.factor(NEO), y = DTP)) + geom_boxplot()
#genes_anno %>% dplyr::filter(DTP != 0) %>% ggplot(aes(x = as.factor(DTP), y = NEO)) + geom_boxplot()

#add PA
for (set_name in names(gs_pa)) {
  genes_anno <- genes_anno %>% mutate(!!set_name := gene %in% gs_pa[[set_name]])
}

#add DEGs
for (set_name in names(gs_degs)) {
  genes_anno <- genes_anno %>% mutate(!!set_name := gene %in% gs_degs[[set_name]])
}

#save
out_path <- paste0(outFolder, "/string.stats.genes.annotation.txt")
write_tsv(x = genes_anno, file = out_path)

#plots

plts <- list()

sig_names <- c("DTP", "NEO")
set_names <- c(names(gs_pa), names(gs_degs))

for (sig_name in sig_names) {
  
  plts[[sig_name]] <- list()
  
  for (set_name in set_names) {
    
    plts[[sig_name]][[set_name]] <- genes_anno %>% 
      ggplot(aes(x = get(set_name), y = get(sig_name))) + 
      geom_boxplot() +
      theme_bw() +
      xlab("") +
      ylab(paste0(sig_name, " (Min. Pervasiveness)")) +
      ggtitle(set_name) +
      theme(plot.title = element_text(size = 8, face = "bold"))
    
  }
  
}

out_F <- paste(outFolder, "/string.stats.genes.annotation.boxplots.pdf", sep = "")
pdf(out_F, width = 2, height = 2.5)
for (sig_name in sig_names) {
  for (set_name in set_names) {
    plot(plts[[sig_name]][[set_name]])
  }
}
dev.off()

genes_anno_summary <- genes_anno %>% 
  pivot_longer(cols = -c("gene", "NEO", "DTP"), names_to = "name", values_to = "in_signature") %>%
  mutate(direction = ifelse(grepl("up", name, ignore.case = TRUE), "up", "down")) %>%
  rowwise() %>%
  mutate(set = gsub("_UP|_DOWN|_top1k-up|_top1k-down", "", name)) %>%
  ungroup()

genes_anno_summary <- genes_anno_summary %>%
  mutate(group = factor(group, levels = rev(c("PA", "First_line", "Resistant", "Second_line", "ER_mutant", "ER_mutant_first_line"))))

#add group
set_groups <- read_tsv("downstream.pervasiveness.string.group.txt")
genes_anno_summary <- genes_anno_summary %>%
  left_join(set_groups, by = "set")

plt_summary_dtp <- genes_anno_summary %>%  
  ggplot(aes(x = group, y = DTP, fill = in_signature)) + 
  geom_boxplot() + 
  facet_wrap(~direction) + 
  coord_flip() + 
  theme_bw() +
  ylab("DTP (Min. Pervasiveness)")

plt_summary_neo <- genes_anno_summary %>%  
  ggplot(aes(x = group, y = NEO, fill = in_signature)) + 
  geom_boxplot() + 
  facet_wrap(~direction) + 
  coord_flip() + 
  theme_bw() +
  ylab("NEO (Min. Pervasiveness)")

out_F <- paste(outFolder, "/string.stats.genes.annotation.boxplots.summary.pdf", sep = "")
pdf(out_F, width = 5, height = 2.5)
plot(plt_summary_dtp)
plot(plt_summary_neo)
dev.off()

