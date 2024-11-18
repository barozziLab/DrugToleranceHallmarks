
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prostate Primary Preparation
#############################

in_folder <- "./prostate_primary/"

#rajan_et_al_2014

in_F <- paste(in_folder, "/rajan_et_al_2014/tax-adt_pre_vs_post.up.txt", sep = "")
in_d <- read_tsv(in_F, col_names = FALSE) %>% pull(1)
m_t2g[["Prostate_primary"]] <- data.frame(gs_name = "rajan_et_al_2014__tax-adt_post_vs_pre_up", gene_symbol = in_d)
in_F <- paste(in_folder, "/rajan_et_al_2014/tax-adt_pre_vs_post.down.txt", sep = "")
in_d <- read_tsv(in_F, col_names = FALSE) %>% pull(1)
m_t2g[["Prostate_primary"]] <- rbind(m_t2g[["Prostate_primary"]], data.frame(gs_name = "rajan_et_al_2014__tax-adt_post_vs_pre_down", gene_symbol = in_d))

#sowalsky_et_al_2018

in_F <- paste(in_folder, "/sowalsky_et_al_2018/limma_results.txt", sep = "")
in_d <- read_tsv(in_F)
t2g_sowalsky <- in_d %>% 
  mutate(gs_name = ifelse(logFC >= log2(1.5) & adj.P.Val <= 0.05, "sowalsky_et_al_2018__adt-intensive_post_vs_ctrl_up", 
                          ifelse(logFC <= -log2(1.5) & adj.P.Val <= 0.05, "sowalsky_et_al_2018__adt-intensive_post_vs_ctrl_down", ""))) %>%
  dplyr::filter(gs_name != "") %>%
  dplyr::select(gs_name, gene_symbol) %>% 
  as.data.frame()
m_t2g[["Prostate_primary"]] <- rbind(m_t2g[["Prostate_primary"]], t2g_sowalsky)

#tewari_et_al_2021

in_F <- paste(in_folder, "/tewari_et_al_2021/deseq2.results.txt", sep = "")
in_d <- read_tsv(in_F)
t2g_tewari <- in_d %>% 
  mutate(gs_name = ifelse(significant_up, "tewari_et_al_2021__adt_resp_vs_non_resp_up",
                          ifelse(significant_down, "tewari_et_al_2021__adt_resp_vs_non_resp_down", ""))) %>%
  dplyr::filter(gs_name != "") %>%
  dplyr::rename(gene_symbol =  gene) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()
m_t2g[["Prostate_primary"]] <- rbind(m_t2g[["Prostate_primary"]], t2g_tewari)

#linder_et_al_2022

in_F <- paste(in_folder, "/linder_et_al_2022/deseq2.results.txt", sep = "")
in_d <- read_tsv(in_F)
t2g_linder <- in_d %>% 
  mutate(gs_name = ifelse(significant_up, "linder_et_al_2022__enz_post_vs_pre_up",
                          ifelse(significant_down, "linder_et_al_2022__enz_post_vs_pre_down", ""))) %>%
  dplyr::filter(gs_name != "") %>%
  dplyr::rename(gene_symbol =  gene) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()
m_t2g[["Prostate_primary"]] <- rbind(m_t2g[["Prostate_primary"]], t2g_linder)

#convert to list (to be queried in GSEA)
pp_list <- list()
for (program_id in unique(m_t2g[["Prostate_primary"]]$gs_name)) {
  pp_list[[program_id]] <- m_t2g[["Prostate_primary"]] %>% dplyr::filter(gs_name == program_id) %>% pull(gene_symbol)
}

#save the full list of markers
out_path <- paste0(out_folder, "/primary.markers.txt")
write_tsv(x = m_t2g[["Prostate_primary"]], file = out_path)


#~~~~~~~~~~~~~~~~~~~~~~~~~
#Pervasiveness Preparation
##########################

perv_sig <- list()

in_folder <- "../literature_meta-analysis/"

in_F <- paste(in_folder, "/DTP-signatures.txt", sep = "")
perv_sig[["DTP"]] <- read_tsv(in_F)
m_t2g[["Breast_pervasinvess"]] <- perv_sig[["DTP"]] %>% 
  dplyr::select(class_2, gene_symbol) %>% 
  dplyr::rename(gs_name = class_2) %>% 
  mutate(gs_name = paste0("DTP_", gs_name)) %>%
  as.data.frame()

in_F <- paste(in_folder, "/NEO-signatures.txt", sep = "")
perv_sig[["NEO"]] <- read_tsv(in_F)
m_t2g_perv_NEO <- perv_sig[["NEO"]] %>% 
  dplyr::select(class_2, gene_symbol) %>% 
  dplyr::rename(gs_name = class_2) %>% 
  mutate(gs_name = paste0("NEO_", gs_name)) %>%
  as.data.frame()
m_t2g[["Breast_pervasinvess"]] <- rbind(m_t2g[["Breast_pervasinvess"]], m_t2g_perv_NEO)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prostate in vitro DEGs Preparation
###################################

# DElegate DEGs

degs_res <- readRDS("../prostate_seurat/deg_results_deseq.rds")

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
m_t2g[["Prostate_invitro_DEGs"]] <- gs_degs_tib


#~~~~~~~~~~~~~~~~~~~~~
#Decipher-seq Programs
######################

m_t2g[["Prostate_invitro_programs"]] <- mrk_tib %>% 
  mutate(program = paste0("Program_", program)) %>% 
  dplyr::rename(gs_name = program) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  as.data.frame()


#~~~~
#GSEA
#####

set_oi <- c("MP_and_PA", "Breast_pervasinvess", "Prostate_invitro_DEGs", "Prostate_invitro_programs")
m_t2g_f <- m_t2g[set_oi]

univ <- do.call(rbind, m_t2g_f) %>% pull(gene_symbol) %>% sort() %>% unique()

res <- list()
res_tib <- c()
for (set_name in names(pp_list)) {
  for (gs_name in names(m_t2g_f)) {
    enricher_res <- enricher(gene = pp_list[[set_name]], TERM2GENE = m_t2g_f[[gs_name]], universe = univ, minGSSize = 10, maxGSSize = 1500)
    if (!is.null(enricher_res)) {
      res[[set_name]] <- enricher_res
      res_tib_i <- res[[set_name]]@result %>% as_tibble() %>% 
        mutate(Gene_set = set_name) %>% 
        mutate(Ontology = gs_name) %>% 
        dplyr::select(Gene_set, Ontology, everything())
      res_tib <- rbind(res_tib, res_tib_i)
    }
  }
}

out_path <- paste0(out_folder, "/primary.enrichr.results.txt")
write_tsv(x = res_tib, file = out_path)

lvls <- c("rajan_et_al_2014__tax-adt_post_vs_pre_up", "rajan_et_al_2014__tax-adt_post_vs_pre_down",
          "linder_et_al_2022__enz_post_vs_pre_up", "linder_et_al_2022__enz_post_vs_pre_down",
          "tewari_et_al_2021__adt_resp_vs_non_resp_up", "tewari_et_al_2021__adt_resp_vs_non_resp_down",
          "sowalsky_et_al_2018__adt-intensive_post_vs_ctrl_up", "sowalsky_et_al_2018__adt-intensive_post_vs_ctrl_down")
res_tib$Gene_set <- factor(res_tib$Gene_set, levels = rev(lvls))

gsea_plts <- list()

for (gs_name in c("Breast_pervasinvess", "MP_and_PA", "Prostate_invitro_DEGs")) {
  
  gsets <- m_t2g_f[[gs_name]] %>% pull(gs_name) %>% unique()
  
  data <- res_tib %>% 
    dplyr::filter(ID %in% gsets) %>% 
    mutate(fraction = Count / top_genes_n) %>%
    dplyr::rename(padjust = p.adjust) %>%
    dplyr::select(Gene_set, ID, Count, fraction, padjust)
  
  padj_cl <- rep("ns", nrow(data))
  padj_cl[data$padjust <= 0.05] <- "adjusted_pval <= 0.05"
  padj_cl[data$padjust <= 0.001] <- "adjusted_pval <= 1e-3"
  padj_cl[data$padjust <= 0.00001] <- "adjusted_pval <= 1e-5"
  padj_cl <- factor(padj_cl, levels = c("adjusted_pval <= 1e-5", "adjusted_pval <= 1e-3", "adjusted_pval <= 0.05", "ns"))
  
  data <- data %>% 
    mutate(padjust_cl = padj_cl) %>%
    mutate(fraction_saturated = ifelse(fraction <= 0.2, fraction, 0.2))
  
  data_wide <- data %>% 
    dplyr::select(Gene_set, ID, fraction) %>% 
    pivot_wider(names_from = ID, values_from = fraction) %>%
    column_to_rownames(var = "Gene_set")
  data_wide[is.na(data_wide)] <- 0
  
  data_wide[is.na(data_wide)] <- 1
  
  #keep only gene sets significant in at least one group
  pval_wide <- data %>% 
    dplyr::select(Gene_set, ID, padjust) %>% 
    pivot_wider(names_from = ID, values_from = padjust) %>%
    column_to_rownames(var = "Gene_set")
  pval_wide[is.na(pval_wide)] <- 1
  gset_which <- apply(pval_wide, 2, min) <= 0.05
  
  data_wide_filt <- data_wide[,gset_which]
  data_filt <- data %>%
    filter(ID %in% colnames(data_wide_filt))
  
  #hc_programs <- hclust(as.dist(1 - cor(t(data_wide_filt), method = "spearman")))
  #fct <- rownames(data_wide_filt)[hc_programs$order]
  #data_filt$Gene_set <- factor(data_filt$Gene_set, levels = fct)
  
  # sort the programs using the 'programs_order' variable
  # from the 'run_factors' step to sort the results
  #data_filt$Gene_set <- factor(data_filt$Gene_set, levels = programs_order)
  
  hc_gsets <- hclust(as.dist(1 - cor(data_wide_filt, method = "spearman")))
  fct <- colnames(data_wide_filt)[hc_gsets$order]
  data_filt$ID <- factor(data_filt$ID, levels = fct)
  
  gsea_plts[[gs_name]] <- ggplot(data_filt, aes(x = ID, y = Gene_set)) + 
    geom_point(aes(size = fraction_saturated, fill = padjust_cl), alpha = 0.75, shape = 21) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_brewer(palette = "Reds", direction = -1) +
    scale_size(range = c(.1, 5), name="Fraction of Genes (Sat.)")
  
}

plt_w <- list()
plt_w[["Breast_pervasinvess"]] <- 6.5
plt_w[["Prostate_invitro_DEGs"]] <- 7
plt_w[["MP_and_PA"]] <- 9

plt_h <- list()
plt_h[["Breast_pervasinvess"]] <- 3.5
plt_h[["Prostate_invitro_DEGs"]] <- 4.5
plt_h[["MP_and_PA"]] <- 3.5

for (gs_name in c("Breast_pervasinvess", "MP_and_PA", "Prostate_invitro_DEGs")) {
  
  out_f <- paste0(out_folder, "/primary.enrichr.results.", gs_name, ".pdf")
  pdf(out_f, height = plt_h[[gs_name]], width = plt_w[[gs_name]])
  plot(gsea_plts[[gs_name]])
  dev.off()
  
}

