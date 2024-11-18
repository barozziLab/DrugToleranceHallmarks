
#~~~~~~~~~~~~~~~~~~
#Factors similarity
###################

col_fun_01 <- colorRamp2(breaks = seq(0, 0.5, 0.05), colors = rev(brewer.pal(11,"Spectral")))

#based on Jaccard index of top genes
programs <- names(mrk_list)
jacc_sim <- c()
for (program_1 in programs) {
  for (program_2 in programs) {
    jacc_sim <- rbind(jacc_sim, cbind(program_1, program_2, jaccard(mrk_list[[program_1]], mrk_list[[program_2]])))
  }
}
jacc_sim <- tibble(program_1 = jacc_sim[,1],
                   program_2 = jacc_sim[,2],
                   JI = as.numeric(jacc_sim[,3]))

jacc_sim_mat <- jacc_sim %>% 
  pivot_wider(names_from = program_1, values_from = JI) %>% 
  column_to_rownames(var = "program_2")

hc <- hclust(as.dist(1 - jacc_sim_mat))

out_path <- paste0(out_folder, "/W.similarity.nGenesPerProgram=", top_genes_n, ".hc.pdf")
pdf(out_path, width = 6, height = 5.5)
plot(hc)
dev.off()

set.seed(123)
hm <- Heatmap(jacc_sim_mat, 
              col = col_fun_01, 
              cluster_rows = TRUE, 
              clustering_distance_rows = "spearman", 
              cluster_columns = TRUE, 
              clustering_distance_columns = "spearman", 
              name = "JI", 
              heatmap_legend_param = list(title = "JI"), 
              row_names_max_width = unit(10, "cm"), 
              border = TRUE)

out_path <- paste0(out_folder, "/W.similarity.nGenesPerProgram=", top_genes_n, ".heatmap.pdf")
pdf(out_path, width = 6, height = 5.5)
draw(hm)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~
#Showing gene vs factors
########################

col_fun <- colorRamp2(breaks = seq(-3.5, -1, 0.25), colors = rev(brewer.pal(11,"Spectral")))
col_fun_scaled <- colorRamp2(c(-2, 0, 2), colors = c("navy", "white", "red"))

## Top Genes

W_topgenes <- W %>% 
  filter(gene %in% mrk_tib$gene_symbol)

W_topgenes_mat <- W_topgenes %>%
  column_to_rownames("gene") %>%
  log10()

hm <- Heatmap(W_topgenes_mat[,programs_order], 
              col = col_fun, 
              cluster_rows = TRUE, 
              clustering_distance_rows = "spearman", 
              cluster_columns = FALSE, 
              name = "W", 
              heatmap_legend_param = list(title = "W"), 
              row_names_max_width = unit(10, "cm"), 
              row_names_gp = gpar(fontsize = 5), 
              border = TRUE)

#note: too big, maybe use only super-top (tried, not that great)

## PA (up)

#exclude ribo

goi <- m_t2g$PA %>% 
  filter(gs_name == "PA_UP") %>% 
  pull(gene_symbol)

#exclude ribo
goi <- goi[!grepl("^RPL|^RPS", goi)]

W_goi <- W %>% 
  filter(gene %in% goi)

W_goi_mat <- W_goi %>%
  column_to_rownames("gene") %>%
  log10()

ha <- HeatmapAnnotation(W = anno_boxplot(W_goi_mat[,programs_order], height = unit(2, "cm")))

set.seed(123)
hm <- Heatmap(W_goi_mat[,programs_order], 
              col = col_fun, 
              cluster_rows = TRUE, 
              cluster_columns = FALSE, 
              row_km = 5,
              name = "W", 
              heatmap_legend_param = list(title = "W"), 
              row_names_max_width = unit(10, "cm"), 
              row_names_gp = gpar(fontsize = 5), 
              border = TRUE,
              bottom_annotation = ha)

out_path <- paste0(out_folder, "/W.PA.up.noribo.heatmap.long.pdf")  
pdf(out_path, width = 6, height = 16)
draw(hm)
dev.off()

set.seed(123)
hm <- Heatmap(W_goi_mat[,programs_order], 
              col = col_fun, 
              cluster_rows = TRUE, 
              cluster_columns = FALSE, 
              row_km = 5,
              name = "W", 
              heatmap_legend_param = list(title = "W"), 
              row_names_max_width = unit(10, "cm"), 
              row_names_gp = gpar(fontsize = 0), 
              column_names_gp = gpar(fontsize = 6), 
              border = TRUE,
              bottom_annotation = ha)

out_path <- paste0(out_folder, "/W.PA.up.noribo.heatmap.short.pdf")  
pdf(out_path, width = 4, height = 6)
draw(hm)
dev.off()

## PA (down)

#exclude ribo

goi <- m_t2g$PA %>% 
  filter(gs_name == "PA_DOWN") %>% 
  pull(gene_symbol)

#exclude ribo
goi <- goi[!grepl("^RPL|^RPS", goi)]

W_goi <- W %>% 
  filter(gene %in% goi)

W_goi_mat <- W_goi %>%
  column_to_rownames("gene") %>%
  log10()

ha <- HeatmapAnnotation(W = anno_boxplot(W_goi_mat[,programs_order], height = unit(2, "cm")))

set.seed(123)
hm <- Heatmap(W_goi_mat[,programs_order], 
              col = col_fun, 
              cluster_rows = TRUE, 
              cluster_columns = FALSE, 
              row_km = 5,
              name = "W", 
              heatmap_legend_param = list(title = "W"), 
              row_names_max_width = unit(10, "cm"), 
              row_names_gp = gpar(fontsize = 5), 
              border = TRUE,
              bottom_annotation = ha)

out_path <- paste0(out_folder, "/W.PA.down.noribo.heatmap.long.pdf")  
pdf(out_path, width = 6, height = 42)
draw(hm)
dev.off()

set.seed(123)
hm <- Heatmap(W_goi_mat[,programs_order], 
              col = col_fun, 
              cluster_rows = TRUE, 
              cluster_columns = FALSE, 
              row_km = 5,
              name = "W", 
              heatmap_legend_param = list(title = "W"), 
              row_names_max_width = unit(10, "cm"), 
              row_names_gp = gpar(fontsize = 0), 
              column_names_gp = gpar(fontsize = 6), 
              border = TRUE,
              bottom_annotation = ha)

out_path <- paste0(out_folder, "/W.PA.down.noribo.heatmap.short.pdf")  
pdf(out_path, width = 4, height = 6)
draw(hm)
dev.off()


#~~~~
#GSEA
#####

univ <- sort(unique(c(colnames(nmf_res$MCF7$W), do.call(rbind, m_t2g_f) %>% pull(gene_symbol))))

res <- list()
res_tib <- c()
for (set_name in names(mrk_list)) {
  for (gs_name in names(m_t2g_f)) {
    enricher_res <- enricher(gene = mrk_list[[set_name]], TERM2GENE = m_t2g_f[[gs_name]], universe = univ, minGSSize = 10, maxGSSize = 1500)
    if (!is.null(enricher_res)) {
      res[[set_name]] <- enricher_res
      res_tib_i <- res[[set_name]]@result %>% as_tibble() %>% 
        mutate(Gene_set = set_name) %>% 
        dplyr::select(Gene_set, everything())
      res_tib <- rbind(res_tib, res_tib_i)
    }
  }
}

out_path <- paste0(out_folder, "/W.enrichr.results.txt") 
write_tsv(x = res_tib, file = out_path)

gsea_plts <- list()
gsea_plts_FL_sorted <- list()
pval_ts <- list()
pval_ts[["MP"]] <- 0.001
pval_ts[["Neo_ET"]] <- 0.05
pval_ts[["Breast_Adult"]] <- 0.05

for (gs_name in c("MP", "Neo_ET", "Breast_Adult")) {

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
  
  #keep only gene sets significant in at least one group
  pval_wide <- data %>% 
    dplyr::select(Gene_set, ID, padjust) %>% 
    pivot_wider(names_from = ID, values_from = padjust) %>%
    column_to_rownames(var = "Gene_set")
  pval_wide[is.na(pval_wide)] <- 1
  gset_which <- apply(pval_wide, 2, min) <= pval_ts[[gs_name]]
  
  data_wide_filt <- data_wide[,gset_which]
  data_filt <- data %>%
    filter(ID %in% colnames(data_wide_filt))
  
  #hc_programs <- hclust(as.dist(1 - cor(t(data_wide_filt), method = "spearman")))
  #fct <- rownames(data_wide_filt)[hc_programs$order]
  #data_filt$Gene_set <- factor(data_filt$Gene_set, levels = fct)
  
  # sort the programs using the 'programs_order' variable
  # from the 'run_factors' step to sort the results
  data_filt$Gene_set <- factor(data_filt$Gene_set, levels = programs_order)
  
  hc_gsets <- hclust(as.dist(1 - cor(data_wide_filt, method = "spearman")))
  fct <- colnames(data_wide_filt)[hc_gsets$order]
  data_filt$ID <- factor(data_filt$ID, levels = fct)
  
  gsea_plts[[gs_name]] <- ggplot(data_filt, aes(x = Gene_set, y = ID)) + 
    geom_point(aes(size = fraction_saturated, fill = padjust_cl), alpha = 0.75, shape = 21) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_brewer(palette = "Reds", direction = -1) +
    scale_size(range = c(.1, 5), name="Fraction of Genes (Sat.)")
  
  gsea_plts_FL_sorted[[gs_name]] <- ggplot(data_filt %>% 
                                             mutate(Gene_set = factor(Gene_set, levels = programs_sorted_FL)), 
                                           aes(x = Gene_set, y = ID)) + 
    geom_point(aes(size = fraction_saturated, fill = padjust_cl), alpha = 0.75, shape = 21) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_brewer(palette = "Reds", direction = -1) +
    scale_size(range = c(.1, 5), name="Fraction of Genes (Sat.)")

}

plt_hs <- list()
plt_hs[["MP"]] <- 4
plt_hs[["Neo_ET"]] <- 1.2
plt_hs[["Breast_Adult"]] <- 2.5

for (gs_name in c("MP", "Neo_ET", "Breast_Adult")) {
  
  out_f <- paste0(out_folder, "/W.enrichr.results.", gs_name, ".nGenesPerProgram=", top_genes_n, ".pdf")
  pdf(out_f, height = plt_hs[[gs_name]], width = 8.5)
  plot(gsea_plts[[gs_name]])
  dev.off()
  
  out_f <- paste0(out_folder, "/W.enrichr.results.", gs_name, ".nGenesPerProgram=", top_genes_n, ".FL-sorted.pdf")
  pdf(out_f, height = plt_hs[[gs_name]], width = 8.5)
  plot(gsea_plts_FL_sorted[[gs_name]])
  dev.off()
  
}


#~~~~~~~~~~
#Statistics
###########

## Fraction of PA up and down fraction per top genes in each program

gsets <- unique(m_t2g[["PA"]]$gs_name)
for (gset in gsets) {
  genes <- m_t2g[["PA"]] %>% dplyr::filter(gs_name == gset) %>% pull(gene_symbol)
  mrk_tib <- mrk_tib %>% mutate(!!gset := gene_symbol %in% genes)
}

out_path <- paste0(out_folder, "/W.PA.genes.txt") 
write_tsv(x = mrk_tib, file = out_path)

mrk_tib_stats <- mrk_tib %>% 
  group_by(program) %>% 
  summarise(n = n(), 
            PA_UP_noRibo_tot = sum(PA_UP_noRibo), 
            PA_UP_noRibo_frac = PA_UP_noRibo_tot / n,
            PA_DOWN_noRibo_tot = sum(PA_DOWN_noRibo), 
            PA_DOWN_noRibo_frac = PA_DOWN_noRibo_tot / n)

mrk_tib_stats <- mrk_tib_stats %>%
  mutate(program = factor(program, levels = programs_order))

plt <- list()

plt[["PA_UP_noRibo"]] <-  ggplot(mrk_tib_stats, aes(x=program, y=PA_UP_noRibo_frac)) + 
  geom_bar(stat = "identity") + 
  theme_bw() +
  ylab("Fraction of genes") +
  xlab("Program") +
  ggtitle("Pre-Adapted up (no ribo genes)") +
  ylim(0,1)

plt[["PA_DOWN_noRibo"]] <-  ggplot(mrk_tib_stats, aes(x=program, y=PA_DOWN_noRibo_frac)) + 
  geom_bar(stat = "identity") + 
  theme_bw() +
  ylab("Fraction of genes") +
  xlab("Program") +
  ggtitle("Pre-Adapted down (no ribo genes)") +
  ylim(0,1)

out_f <- paste0(out_folder, "/W.PA.barplot.pdf")
pdf(out_f, height = 5, width = 7)
grid.arrange(plt[["PA_UP_noRibo"]], plt[["PA_DOWN_noRibo"]], nrow = 2)
dev.off()

#sorted by programs_sorted_FL

plt <- list()

plt[["PA_UP_noRibo"]] <-  ggplot(mrk_tib_stats %>% 
                                   mutate(program = factor(program, levels = programs_sorted_FL)), 
                                 aes(x=program, y=PA_UP_noRibo_frac)) + 
  geom_bar(stat = "identity") + 
  theme_bw() +
  ylab("Fraction of genes") +
  xlab("Program") +
  ggtitle("Pre-Adapted up (no ribo genes)") +
  ylim(0,1)

plt[["PA_DOWN_noRibo"]] <-  ggplot(mrk_tib_stats %>%
                                     mutate(program = factor(program, levels = programs_sorted_FL)),
                                   aes(x=program, y=PA_DOWN_noRibo_frac)) + 
  geom_bar(stat = "identity") + 
  theme_bw() +
  ylab("Fraction of genes") +
  xlab("Program") +
  ggtitle("Pre-Adapted down (no ribo genes)") +
  ylim(0,1)

out_f <- paste0(out_folder, "/W.PA.barplot.FL-sorted.pdf")
pdf(out_f, height = 5, width = 7)
grid.arrange(plt[["PA_UP_noRibo"]], plt[["PA_DOWN_noRibo"]], nrow = 2)
dev.off()


#~~~~~
#FGSEA
######

#proper GSEA score using ranked genes (fgsea package), so that no need to threshold on the coefficients

gs_name <- "MP_and_PA"

sets_oi <- m_t2g[[gs_name]]
sets_oi_list <- list()
for (set in unique(sets_oi$gs_name)) {
  goi <- sets_oi %>% filter(gs_name == set) %>% pull(gene_symbol)
  sets_oi_list[[set]] <- goi
}

program_ids <- W %>% colnames()
program_ids <- program_ids[program_ids != "gene"]

fgseaRes <- list()

for (program_id in program_ids) {
  
  neo_ranked_genes_lst <- W %>% arrange(!!program_id) %>% pull(!!program_id)
  names(neo_ranked_genes_lst) <- W %>% arrange(!!program_id) %>% pull(gene)
  neo_ranked_genes_lst <- log10(neo_ranked_genes_lst)
  neo_ranked_genes_lst <- sort(neo_ranked_genes_lst)
  
  #excluding the "missing" ones
  neo_ranked_genes_lst <- neo_ranked_genes_lst[neo_ranked_genes_lst >= -15]
  
  fgseaRes[[program_id]] <- fgsea(pathways = sets_oi_list,
                                  stats    = neo_ranked_genes_lst,
                                  eps      = 0.0,
                                  minSize  = 10,
                                  maxSize  = 1500,
                                  nPermSimple = 100000)
  
}

#merge
fgseaRes_long <- c()
for (program_id in names(fgseaRes)) {
  fgseaRes_long <- rbind(fgseaRes_long, cbind(program_id, fgseaRes[[program_id]]))
}
fgseaRes_long_full <- fgseaRes_long
fgseaRes_long <- fgseaRes_long_full %>% select(-leadingEdge)

#save resuts
out_path <- paste0(out_folder, "/W.fgsea.results.rds") 
saveRDS(object = fgseaRes, file = out_path)
out_path <- paste0(out_folder, "/W.fgsea.results.txt") 
write_tsv(x = fgseaRes_long, file = out_path)

fgseaRes_wide <- fgseaRes_long %>% 
  dplyr::select(pathway, program_id, NES) %>% 
  pivot_wider(names_from = program_id, values_from = NES) %>%
  column_to_rownames(var = "pathway")
fgseaRes_wide[is.na(fgseaRes_wide)] <- 0

zeros <- apply(fgseaRes_wide == 0, 1, sum) / ncol(fgseaRes_wide)
fgseaRes_wide <- fgseaRes_wide[!zeros > 0.9, ]

col_fun_scaled_new <- colorRamp2(c(-4, 0, 4), colors = c("navy", "white", "red"))

hm <- Heatmap(fgseaRes_wide[,programs_order], 
              col = col_fun_scaled_new, 
              cluster_rows = TRUE, 
              clustering_distance_rows = "spearman", 
              cluster_columns = FALSE, 
              name = "NES", 
              heatmap_legend_param = list(title = "NES"), 
              row_names_max_width = unit(10, "cm"), 
              border = TRUE)

out_f <- paste0(out_folder, "/W.fgsea.results.", gs_name, ".pdf")
pdf(out_f, height = 6, width = 9.5)
draw(hm)
dev.off()

#also sorted on the FL subtracted values

hm <- Heatmap(fgseaRes_wide[,programs_sorted_FL], 
              col = col_fun_scaled_new, 
              cluster_rows = TRUE, 
              clustering_distance_rows = "spearman", 
              cluster_columns = FALSE, 
              name = "NES", 
              heatmap_legend_param = list(title = "NES"), 
              row_names_max_width = unit(10, "cm"), 
              border = TRUE)

out_f <- paste0(out_folder, "/W.fgsea.results.", gs_name, ".FL-sorted.pdf")
pdf(out_f, height = 6, width = 9.5)
draw(hm)
dev.off()

#PA signatures barplots

nes_sat <- 5

bplots <- list()

sigs <- c("PA_UP", "PA_UP_noRibo", "PA_DOWN", "PA_DOWN_noRibo")

#original ordering

for (sig in sigs) {
  
  bplots[[sig]] <- fgseaRes_long %>% 
    dplyr::filter(pathway == sig) %>% 
    mutate(program_id = factor(program_id, levels = rev(programs_order))) %>% 
    mutate(NES = ifelse(NES >= nes_sat, nes_sat, ifelse(NES <= -nes_sat, -nes_sat, NES))) %>%
    mutate(NES = ifelse(is.na(NES), 0, NES)) %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    ggplot(aes(x = program_id, y = NES, fill = -log10(padj))) + 
    geom_bar(stat="identity") + 
    coord_flip() +
    theme_bw() + 
    ylim(-nes_sat, nes_sat) +
    ggtitle(sig)
  
}

out_f <- paste0(out_folder, "/W.fgsea.results.PA.barplot.pdf")
pdf(out_f, height = 4.5, width = 6)
plot(bplots[["PA_UP"]] + bplots[["PA_DOWN"]])
plot(bplots[["PA_UP_noRibo"]] + bplots[["PA_DOWN_noRibo"]])
dev.off()

#FL ordering

bplots <- list()

for (sig in sigs) {
  
  bplots[[sig]] <- fgseaRes_long %>% 
    dplyr::filter(pathway == sig) %>% 
    mutate(program_id = factor(program_id, levels = rev(programs_sorted_FL))) %>% 
    mutate(NES = ifelse(NES >= nes_sat, nes_sat, ifelse(NES <= -nes_sat, -nes_sat, NES))) %>%
    mutate(NES = ifelse(is.na(NES), 0, NES)) %>%
    mutate(padj = ifelse(is.na(padj), 1, padj)) %>%
    ggplot(aes(x = program_id, y = NES, fill = -log10(padj))) + 
    geom_bar(stat="identity") + 
    coord_flip() +
    theme_bw() + 
    ylim(-nes_sat, nes_sat) +
    ggtitle(sig)
  
}

out_f <- paste0(out_folder, "/W.fgsea.results.PA.barplot.FL-sorted.pdf")
pdf(out_f, height = 4.5, width = 6)
plot(bplots[["PA_UP"]] + bplots[["PA_DOWN"]])
plot(bplots[["PA_UP_noRibo"]] + bplots[["PA_DOWN_noRibo"]])
dev.off()

#PA UP no ribo genes

goi <- c()
for (i in 1:nrow(fgseaRes_long_full)) {
  flag <- as.vector(fgseaRes_long_full[i,"pathway"] == "PA_UP_noRibo")
  if (flag) {
    res_flag <- fgseaRes_long_full[i,c("program_id", "leadingEdge")]
    goi <- rbind(goi, cbind(res_flag$program_id, res_flag$leadingEdge[[1]]))
  }
}
genes_paupnoribo <- tibble(program = goi[,1], gene = goi[,2])
genes_paupnoribo$seen <- 1
genes_paupnoribo <- genes_paupnoribo %>% pivot_wider(names_from = program, values_from = seen)
genes_paupnoribo[is.na(genes_paupnoribo)] <- 0

