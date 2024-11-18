
source("downstream.libs.R")

outFolder <- "downstream.pervasiveness.results/"


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
t2g[["Breast_DEGs"]] <- gs_degs_tib

## Msigdb hallmarks

msig_extract <- msigdbr(species = "Homo sapiens", category = "H")
t2g[["MSigDB_Hallmark"]] <- msig_extract %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

## Msigdb reactome

msig_extract <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
t2g[["MSigDB_Reactome"]] <- msig_extract %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

## Msigdb KEGG (includes signaling pathways)

msig_extract <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
t2g[["MSigDB_KEGG"]] <- msig_extract %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

## Dorothea for regulons

dorothea_hs_f <- dorothea_hs %>% dplyr::filter(confidence != "E")
t2g[["Dorothea"]] <- data.frame(gs_name = dorothea_hs_f %>% pull(tf),
                                gene_symbol = dorothea_hs_f %>% pull(target),
                                stringsAsFactors = FALSE)

## Meta-programs Tirosh

MP <- readRDS("../data/MP.rds")
t2g[["MP"]] <- MP$hs %>% dplyr::rename(gene_symbol = gene) %>% as.data.frame()


#~~~~~~~~
#Analysis
#########

## fGSEA

fgseaRes <- list()

for (gs_name in names(t2g)) {
  
  sets_oi <- t2g[[gs_name]]
  sets_oi_list <- list()
  for (set in unique(sets_oi$gs_name)) {
    goi <- sets_oi %>% filter(gs_name == set) %>% pull(gene_symbol)
    sets_oi_list[[set]] <- goi
  }
  
  program_ids <- names(perv_sig)
  
  fgseaRes[[gs_name]] <- list()
  
  for (program_id in program_ids) {
    
    ranked_genes_lst <- perv_sig[[program_id]] %>% pull(n_sig)
    names(ranked_genes_lst) <- perv_sig[[program_id]] %>% pull(gene_symbol)
    ranked_genes_lst <- sort(ranked_genes_lst, decreasing = TRUE)
    
    fgseaRes[[gs_name]][[program_id]] <- fgsea(pathways = sets_oi_list,
                                               stats    = ranked_genes_lst,
                                               eps      = 0.0,
                                               minSize  = 10,
                                               maxSize  = 1500,
                                               nPermSimple = 10000)
    
  }
  
}

#PA fgsea plots

gs_name <- "Pre-Adapted"

sets_oi <- t2g[[gs_name]]
sets_oi_list <- list()
for (set in unique(sets_oi$gs_name)) {
  goi <- sets_oi %>% filter(gs_name == set) %>% pull(gene_symbol)
  sets_oi_list[[set]] <- goi
}

out_path <- paste0(outFolder, "/fgsea.results.PA.pdf")
pdf(out_path, width = 5, height = 3)

for (program_id in c("NEO", "DTP")) {
  
  sets <- fgseaRes[[gs_name]][[program_id]]$pathway
  
  for (set in sets) {
  
    NES <- fgseaRes[[gs_name]][[program_id]] %>% filter(pathway == set) %>% pull(NES) %>% round(2)
    padj <- fgseaRes[[gs_name]][[program_id]] %>% filter(pathway == set) %>% pull(padj) %>% formatC(2)
    plot_title <- paste(program_id, "; ", set, "; NES = ", NES, "; adjusted p-value = ", padj, sep = "")
    res <- plotEnrichment(sets_oi_list[[set]], ranked_genes_lst) + 
      labs(title=plot_title) +
      theme(plot.title = element_text(size=8))
    plot(res)
    
  }

}

dev.off()

#merge
fgseaRes_long <- c()
for (gs_name in names(fgseaRes)) {
  for (program_id in names(fgseaRes[[gs_name]])) {
    fgseaRes_long <- rbind(fgseaRes_long, cbind(program_id, gs_name, fgseaRes[[gs_name]][[program_id]]))
  }
}
fgseaRes_long_wLeadEdges <- fgseaRes_long %>% as_tibble()
fgseaRes_long <- fgseaRes_long_wLeadEdges %>% select(-leadingEdge)

#save results
out_path <- paste0(outFolder, "/fgsea.results.txt") 
write_tsv(x = fgseaRes_long_wLeadEdges %>% rowwise() %>% mutate(leadingEdge = paste0(sort(leadingEdge), collapse = ",")), file = out_path)

#DTP vs NEO enrichment

fgsea_dtp <- fgseaRes_long %>% dplyr::filter(program_id == "DTP") %>% select(-c("program_id"))
colnames(fgsea_dtp) <- paste0("DTP_", colnames(fgsea_dtp))
fgsea_dtp <- fgsea_dtp %>% dplyr::rename(gs_name = DTP_gs_name, pathway = DTP_pathway)

fgsea_neo <- fgseaRes_long %>% dplyr::filter(program_id == "NEO") %>% select(-c("program_id"))
colnames(fgsea_neo) <- paste0("NEO_", colnames(fgsea_neo))
fgsea_neo <- fgsea_neo %>% dplyr::rename(gs_name = NEO_gs_name, pathway = NEO_pathway)

fgsea_merge <- fgsea_neo %>% full_join(fgsea_dtp, by = c("gs_name", "pathway"))
fgsea_merge$DTP_NES[is.na(fgsea_merge$DTP_NES)] <- 0
fgsea_merge$DTP_padj[is.na(fgsea_merge$DTP_padj)] <- 1
fgsea_merge$NEO_NES[is.na(fgsea_merge$NEO_NES)] <- 0
fgsea_merge$NEO_padj[is.na(fgsea_merge$NEO_padj)] <- 1

#save merge results
out_path <- paste0(outFolder, "/fgsea.results.merge.txt") 
write_tsv(x = fgsea_merge, file = out_path)

## Summary

#plot summary of the results of fgsea on Breast_DEGs

d_input <- fgseaRes_long %>% 
  dplyr::filter(gs_name == "Breast_DEGs") 

pathways_order <- d_input %>% 
  dplyr::select(program_id, pathway, NES) %>% 
  pivot_wider(names_from = program_id, values_from = NES) %>% 
  mutate(strength = DTP + NEO) %>% 
  arrange(-strength) %>% 
  pull(pathway)

d_input <- d_input %>%
  mutate(pathway = factor(pathway, levels = pathways_order))

plt <- d_input %>%
  ggplot(aes(fill=program_id, y=NES, x=pathway)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() + 
  theme_bw()

out_F <- paste(outFolder, "/fgsea.results.comparison.Breast_DEGs.pdf", sep = "")
pdf(out_F, width = 7.5, height = 9)
plot(plt)
dev.off()

#plot summary of the results of fgsea on Dorothea (significant only)

sig_t <- 0.005
genen_t <- 20

d_input <- fgseaRes_long %>% 
  dplyr::filter(gs_name == "Dorothea" & padj <= sig_t & size >= genen_t) 

pathways_order <- d_input %>% 
  dplyr::select(program_id, pathway, NES) %>% 
  pivot_wider(names_from = program_id, values_from = NES) %>% 
  mutate(strength = DTP + NEO) %>% 
  arrange(-strength) %>% 
  pull(pathway)

d_input <- d_input %>%
  mutate(pathway = factor(pathway, levels = pathways_order))

plt <- d_input %>%
  ggplot(aes(fill=program_id, y=NES, x=pathway)) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() + 
  theme_bw()

out_F <- paste(outFolder, "/fgsea.results.comparison.Dorothea.pdf", sep = "")
pdf(out_F, width = 5, height = 10)
plot(plt)
dev.off()

## Network of the most enriched terms

#based on similarity (using jaccard index) of leading edge genes

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

ji_calc <- function(res_lst) {
  ji_mat <- c()
  names_oi <- names(res_lst) 
  for (name_oi_1 in names_oi) {
    for (name_oi_2 in names_oi) {
      ji_mat <- rbind(ji_mat, c(name_oi_1, name_oi_2, jaccard(res_lst[[name_oi_1]], res_lst[[name_oi_2]])))
    }
  }
  ji_tib <- tibble(id_1 = ji_mat[,1],
                   id_2 = ji_mat[,2],
                   JI = as.numeric(ji_mat[,3]))
  ji_df <- ji_tib %>% pivot_wider(names_from = id_2, values_from = JI) %>% column_to_rownames(var = "id_1")
  return(list(ji_tib, ji_df))
}

sig_t <- 0.005
genen_t <- 20

ji_res <- list()
ji_anno <- list()

#NEO, all

ids <- fgsea_merge %>% dplyr::filter(NEO_padj <= sig_t & NEO_size >= genen_t) %>% pull(pathway)

res <- fgseaRes_long_wLeadEdges %>% dplyr::filter(pathway %in% ids & program_id == "NEO")
res_lst <- res %>% pull(leadingEdge)
names(res_lst) <- res %>% pull(pathway)

ji_res[["NEO-all"]] <- ji_calc(res_lst)

#NEO, regulators

ids <- fgsea_merge %>% dplyr::filter(NEO_padj <= sig_t & NEO_size >= genen_t) %>% dplyr::filter(gs_name == "Dorothea") %>% pull(pathway)

res <- fgseaRes_long_wLeadEdges %>% dplyr::filter(pathway %in% ids & program_id == "NEO")
res_lst <- res %>% pull(leadingEdge)
names(res_lst) <- res %>% pull(pathway)

ji_res[["NEO-regulators"]] <- ji_calc(res_lst)

#DTP, all

ids <- fgsea_merge %>% dplyr::filter(DTP_padj <= sig_t & DTP_size >= genen_t) %>% pull(pathway)

res <- fgseaRes_long_wLeadEdges %>% dplyr::filter(pathway %in% ids & program_id == "DTP")
res_lst <- res %>% pull(leadingEdge)
names(res_lst) <- res %>% pull(pathway)

ji_res[["DTP-all"]] <- ji_calc(res_lst)

#DTP, regulators

ids <- fgsea_merge %>% dplyr::filter(DTP_padj <= sig_t & DTP_size >= genen_t) %>% dplyr::filter(gs_name == "Dorothea") %>% pull(pathway)

res <- fgseaRes_long_wLeadEdges %>% dplyr::filter(pathway %in% ids & program_id == "DTP")
res_lst <- res %>% pull(leadingEdge)
names(res_lst) <- res %>% pull(pathway)

ji_res[["DTP-regulators"]] <- ji_calc(res_lst)

#NEO or DTP, regulators

ids_neo <- fgsea_neo %>% 
  dplyr::filter(NEO_padj <= sig_t & NEO_size >= genen_t) %>% 
  dplyr::filter(gs_name == "Dorothea") %>% pull(pathway)

ids_dtp <- fgsea_dtp %>% 
  dplyr::filter(DTP_padj <= sig_t & DTP_size >= genen_t) %>% 
  dplyr::filter(gs_name == "Dorothea") %>% pull(pathway)

res <- fgseaRes_long_wLeadEdges %>% dplyr::filter((pathway %in% ids_neo & program_id == "NEO") | (pathway %in% ids_dtp & program_id == "DTP"))
res_lst <- res %>% pull(leadingEdge)
names(res_lst) <- res %>% rowwise() %>% mutate(id = paste0(pathway, "_", program_id)) %>% pull(id)

ji_res[["Both-regulators"]] <- ji_calc(res_lst)
ji_anno[["Both-regulators"]] <- res %>% pull(NES)

#Heat maps

set.seed(123)
HM <- Heatmap(ji_res[["NEO-all"]][[2]],
              heatmap_legend_param = list(title = "JI"),
              col = colorRamp2(c(0,0.2), c("white", "red")),
              row_names_gp = grid::gpar(fontsize = 7),
              column_names_gp = grid::gpar(fontsize = 7),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(8, "cm"),
              column_dend_height = unit(8, "cm"),
              row_names_max_width = unit(18, "cm"),
              column_names_max_height = unit(18, "cm"),
              border = TRUE,
              clustering_distance_rows = "pearson",
              clustering_distance_columns = "pearson")

out_F <- paste(outFolder, "/fgsea.results.NEO_top_all.JI.pdf", sep = "")
pdf(out_F, width = 50, height = 50)
draw(HM)
dev.off()

set.seed(123)
HM <- Heatmap(ji_res[["NEO-regulators"]][[2]],
              heatmap_legend_param = list(title = "JI"),
              col = colorRamp2(c(0,0.2), c("white", "red")),
              row_names_gp = grid::gpar(fontsize = 7),
              column_names_gp = grid::gpar(fontsize = 7),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(1, "cm"),
              column_dend_height = unit(1, "cm"),
              row_names_max_width = unit(2, "cm"),
              column_names_max_height = unit(2, "cm"),
              border = TRUE,
              clustering_distance_rows = "pearson",
              clustering_distance_columns = "pearson")

out_F <- paste(outFolder, "/fgsea.results.NEO_top_regulators.JI.pdf", sep = "")
pdf(out_F, width = 6, height = 5)
draw(HM)
dev.off()

set.seed(123)
HM <- Heatmap(ji_res[["DTP-all"]][[2]],
              heatmap_legend_param = list(title = "JI"),
              col = colorRamp2(c(0,0.2), c("white", "red")),
              row_names_gp = grid::gpar(fontsize = 7),
              column_names_gp = grid::gpar(fontsize = 7),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(8, "cm"),
              column_dend_height = unit(8, "cm"),
              row_names_max_width = unit(18, "cm"),
              column_names_max_height = unit(18, "cm"),
              border = TRUE,
              clustering_distance_rows = "pearson",
              clustering_distance_columns = "pearson")

out_F <- paste(outFolder, "/fgsea.results.DTP_top_all.JI.pdf", sep = "")
pdf(out_F, width = 50, height = 50)
draw(HM)
dev.off()

set.seed(123)
HM <- Heatmap(ji_res[["DTP-regulators"]][[2]],
              heatmap_legend_param = list(title = "JI"),
              col = colorRamp2(c(0,0.2), c("white", "red")),
              row_names_gp = grid::gpar(fontsize = 7),
              column_names_gp = grid::gpar(fontsize = 7),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(1, "cm"),
              column_dend_height = unit(1, "cm"),
              row_names_max_width = unit(2, "cm"),
              column_names_max_height = unit(2, "cm"),
              border = TRUE,
              clustering_distance_rows = "pearson",
              clustering_distance_columns = "pearson")

out_F <- paste(outFolder, "/fgsea.results.DTP_top_regulators.JI.pdf", sep = "")
pdf(out_F, width = 9, height = 8)
draw(HM)
dev.off()

ha_cols <- list(NES =  colorRamp2(c(-3,0,3), c("darkgreen", "white", "purple")))

ha <- HeatmapAnnotation(NES = ji_anno[["Both-regulators"]], col = ha_cols)
row_ha <- rowAnnotation(NES = ji_anno[["Both-regulators"]], col = ha_cols)

set.seed(123)
HM <- Heatmap(ji_res[["Both-regulators"]][[2]],
              heatmap_legend_param = list(title = "JI"),
              col = colorRamp2(c(0,0.2), c("white", "red")),
              row_names_gp = grid::gpar(fontsize = 7),
              column_names_gp = grid::gpar(fontsize = 7),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(1, "cm"),
              column_dend_height = unit(1, "cm"),
              row_names_max_width = unit(2, "cm"),
              column_names_max_height = unit(2, "cm"),
              border = TRUE,
              top_annotation = ha,
              right_annotation = row_ha,
              clustering_distance_rows = "pearson",
              clustering_distance_columns = "pearson")

out_F <- paste(outFolder, "/fgsea.results.both_top_regulators.JI.pdf", sep = "")
pdf(out_F, width = 13, height = 12)
draw(HM)
dev.off()

# Bubble plot of the top enriched

top_n <-  20

#add p-value brackets
data <- fgseaRes_long
padj_cl <- rep("ns", nrow(data))
padj_cl[data$padj <= 0.05] <- "adjusted_pval <= 0.05"
padj_cl[data$padj <= 0.001] <- "adjusted_pval <= 1e-3"
padj_cl[data$padj <= 1e-5] <- "adjusted_pval <= 1e-5"
padj_cl <- factor(padj_cl, levels = c("ns", "adjusted_pval <= 0.05", "adjusted_pval <= 1e-3", "adjusted_pval <= 1e-5"))
data <- data %>% mutate(padjust_cl = padj_cl)

groups_oi <- data$gs_name %>% unique()

plts <- list()

for (group_oi in groups_oi) {
  
  data_oi <- data %>% dplyr::filter(gs_name == group_oi)
  data_oi %>% dplyr::filter(program_id == "DTP") %>% arrange(-abs(NES))
  top_dtp <- data_oi %>% dplyr::filter(program_id == "DTP") %>% arrange(-abs(NES)) %>% pull(pathway) %>% head(top_n)
  top_neo <- data_oi %>% dplyr::filter(program_id == "NEO") %>% arrange(-abs(NES)) %>% pull(pathway) %>% head(top_n)
  data_oi <- data_oi %>% dplyr::filter(pathway %in% c(top_dtp, top_neo))
  
  #sort by strongest combined score
  pathways_sorted <- data_oi %>% mutate(score = NES * -log10(padj)) %>% arrange(-score) %>% pull(pathway) %>% unique()
  data_oi$pathway <- factor(data_oi$pathway, levels = pathways_sorted)
  
  plts[[group_oi]] <- ggplot(data_oi, aes(x = program_id, y = pathway)) +
    geom_point(aes(size = padjust_cl, fill = NES), alpha = 0.75, shape = 21) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_distiller(palette = "Spectral")

}

#Add merged bubble plot Pre-Adapted + MSigDB_Hallmark
data_oi <- data %>% dplyr::filter(gs_name %in% c("Pre-Adapted", "MSigDB_Hallmark"))
data_oi %>% dplyr::filter(program_id == "DTP") %>% arrange(-abs(NES))
top_dtp <- data_oi %>% dplyr::filter(program_id == "DTP") %>% arrange(-abs(NES)) %>% pull(pathway) %>% head(top_n)
top_neo <- data_oi %>% dplyr::filter(program_id == "NEO") %>% arrange(-abs(NES)) %>% pull(pathway) %>% head(top_n)
data_oi <- data_oi %>% dplyr::filter(pathway %in% c(top_dtp, top_neo))

#sort by strongest combined score
pathways_sorted <- data_oi %>% mutate(score = NES * -log10(padj)) %>% arrange(-score) %>% pull(pathway) %>% unique()
data_oi$pathway <- factor(data_oi$pathway, levels = pathways_sorted)

plts[["Pre-Adapted__MSigDB_Hallmark"]] <- ggplot(data_oi, aes(x = program_id, y = pathway)) +
  geom_point(aes(size = padjust_cl, fill = NES), alpha = 0.75, shape = 21) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_distiller(palette = "Spectral")

plts_width <- list()
plts_width[["Pre-Adapted"]] <- 4
plts_width[["Breast_DEGs"]] <- 6.5
plts_width[["MSigDB_Hallmark"]] <- 6.5
plts_width[["MSigDB_Reactome"]] <- 8.5
plts_width[["MSigDB_KEGG"]] <- 7
plts_width[["Dorothea"]] <- 3.5
plts_width[["MP"]] <- 4.5
plts_width[["Pre-Adapted__MSigDB_Hallmark"]] <- 6.5

plts_height <- list()
plts_height[["Pre-Adapted"]] <- 4
plts_height[["Breast_DEGs"]] <- 7
plts_height[["MSigDB_Hallmark"]] <- 6
plts_height[["MSigDB_Reactome"]] <- 8
plts_height[["MSigDB_KEGG"]] <- 7
plts_height[["Dorothea"]] <- 7
plts_height[["MP"]] <- 6
plts_height[["Pre-Adapted__MSigDB_Hallmark"]] <- 6.5

for (group_oi in c(groups_oi, "Pre-Adapted__MSigDB_Hallmark")) {
  
  out_F <- paste(outFolder, "/fgsea.results.comparison.bubbleplot.", group_oi, ".pdf", sep = "")
  pdf(out_F, width = plts_width[[group_oi]], height = plts_height[[group_oi]])
  plot(plts[[group_oi]])
  dev.off()
  
}

