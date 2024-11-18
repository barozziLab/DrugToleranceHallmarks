
out_folder <- "decoupleR.stats.run_subtracted-paired/"
cmd <- paste0("mkdir ", out_folder)
system(cmd)


#~~~~~~~~~~~~~~~~
#Data preparation
#################

# Overall, annotated dataset (TF by sample)

d_sub_wide <- d_sub %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  left_join(d_neo, by = "gene") %>%
  dplyr::rename(neoadj_score = score, 
                neoadj_pvalue = p_value, 
                neoadj_pvalue_BH = p_value_BH,
                neoadj_lasso_est = lasso_est,
                neoadj_rf_imp = rf_imp) %>%
  select(-c(condition)) %>%
  mutate(neoadj_class = ifelse(neoadj_score > 0 & neoadj_pvalue_BH <= 0.05, "up", ifelse(neoadj_score < 0 & neoadj_pvalue_BH <= 0.05, "down", "nc")))

#add average first and second line
d_sub_wide <- d_sub_wide %>% 
  rowwise() %>% 
  mutate(first_line = mean(c(`WM-2d`, 
                             `Tam-2d`, 
                             `Fulv-2d`))) %>%
  mutate(second_line = mean(c(`Lted_Tam-7d`, 
                              `Lted_Fulv-7d`, 
                              `Lted_Palbo-7d`, 
                              `Lted_CDK7i-7d`, 
                              `TamR_WM-7d`, 
                              `TamR_Fulv-7d`, 
                              `TamR_Palbo-7d`, 
                              `TamR_CDK7i-7d`)))

# Analyses based on differences in any of the conditions

cols_w <- d_sub_long$name %>% unique() %>% as.character()

res <- list()

res_up <- c()
res_down <- c()

for (col_w in cols_w) {
  
  res[[col_w]] <- data_filter(d_sub_long = d_sub_long, 
                              tf_expr_frac = tf_expr_frac, 
                              sample_w = col_w, 
                              filt_tf_expr = filt_tf_expr,
                              filt_tf_act = filt_tf_act)
  
  reg_up <- res[[col_w]] %>% dplyr::filter(group == "up") %>% pull(gene) %>% unique() %>% as.character()
  res_up[[col_w]] <- c(res_up[[col_w]], reg_up) %>% unique()
  
  reg_down <- res[[col_w]] %>% dplyr::filter(group == "down") %>% pull(gene) %>% unique() %>% as.character()
  res_down[[col_w]] <- c(res_down[[col_w]], reg_down) %>% unique()
  
}

out_file <- paste0(out_folder, "/", "data.all.res.rds")
saveRDS(object = res, file = out_file)


#~~~~~~~~~~~~~~~~
#Overall analysis
#################

#prepare lists of sample ids, per group
ids <- list()
#first line
ids[["First_line"]] <- metadata_ordered %>% dplyr::filter(class == "first_line") %>% pull(cell_type_2) %>% unique()
#mutant
ids[["ESR1_mutant"]] <- "D538G_RM"
#second line, no CDK7i
ids[["Second_line"]] <- metadata_ordered %>% dplyr::filter(class == "second_line" & grepl("7d", cell_type_2)) %>% pull(cell_type_2) %>% unique()
w <- grepl("CDK7i", ids[["Second_line"]])
ids[["Second_line"]] <- ids[["Second_line"]][!w]
#resistant
ids[["Resistant"]] <- metadata_ordered %>% dplyr::filter(class == "resistant") %>% pull(cell_type_2) %>% unique()

ids_all <- unlist(ids)  %>% as.character() %>% sort() %>% unique()

#single table, with a column specifying grouping

d_goi <- d_sub_wide %>% 
  select(gene, all_of(ids_all)) %>% 
  pivot_longer(cols = -c("gene"), names_to = "sample")

sample_class <- rep(NA, nrow(d_goi))
for (id in names(ids)) {
  sample_class[d_goi$sample %in% ids[[id]]] <- id
}
d_goi <- d_goi %>%
  mutate(sample_cls = sample_class)

#subset on specific sets of ids (any significant, at the moment)

#any up-regulation
all_up <- unlist(res_up %>% unlist()) %>% as.character() %>% sort() %>% unique()
#any down-regulation
all_down <- unlist(res_down %>% unlist()) %>% as.character() %>% sort() %>% unique()

d_goi_summary <- d_goi %>% 
  dplyr::filter(gene %in% all_up | gene %in% all_down) %>% 
  group_by(gene, sample_cls) %>% summarise(mean_value = mean(value)) %>%
  ungroup()
d_goi_summary_w <- d_goi_summary %>% 
  pivot_wider(names_from = sample_cls, values_from = mean_value)
d_goi_summary_w_verb <- d_goi_summary_w %>%
  left_join(d_sub_wide, by = "gene") %>%
  select(-c(first_line, second_line))

out_file <- paste0(out_folder, "/", "data.all.wide.anno.txt")
write_tsv(x = d_goi_summary_w_verb, file = out_file)

#Same file but listing all the TFs (for consolidation purposes)

d_goi_summary_complete <- d_goi %>% 
  group_by(gene, sample_cls) %>% summarise(mean_value = mean(value)) %>%
  ungroup() %>% 
  pivot_wider(names_from = sample_cls, values_from = mean_value) %>%
  left_join(d_sub_wide, by = "gene") %>%
  select(-c(first_line, second_line))

out_file <- paste0(out_folder, "/", "data.all.wide.anno.complete_info.txt")
write_tsv(x = d_goi_summary_complete, file = out_file)

#Bar plots of top scoring

n_top <- 20

#Sorted on ESR1 mutant, First_line, etc

for (w_col in c("ESR1_mutant", "First_line", "Second_line", "Resistant")) {
  
  goi <- d_goi_summary_w_verb %>% arrange(-get(w_col)) %>% pull(gene) %>% head(n_top) %>% rev()
  
  plt_up <- d_goi_summary %>% 
    mutate(mean_value = ifelse(mean_value >= 5, 5, ifelse(mean_value <= -3, 3, mean_value))) %>%
    dplyr::filter(gene %in% goi) %>% 
    mutate(gene = factor(gene, levels = goi)) %>%
    ggplot(aes(x=gene, y=mean_value)) +
    geom_bar(stat="identity") + 
    ylim(-3, 5) + 
    xlab("Top Regulator (up)") +
    ylab("Mean Activity (vs baseline)") +
    facet_wrap(~sample_cls, ncol = 4) +
    coord_flip() + 
    theme_bw()
  
  goi <- d_goi_summary_w_verb %>% arrange(get(w_col)) %>% pull(gene) %>% head(n_top)
  
  plt_down <- d_goi_summary %>% 
    mutate(mean_value = ifelse(mean_value >= 5, 5, ifelse(mean_value <= -3, 3, mean_value))) %>%
    dplyr::filter(gene %in% goi) %>% 
    mutate(gene = factor(gene, levels = goi)) %>%
    ggplot(aes(x=gene, y=mean_value)) +
    geom_bar(stat="identity") + 
    ylim(-3, 5) + 
    xlab("Top Regulator (down)") +
    ylab("Mean Activity (vs baseline)") +
    facet_wrap(~sample_cls, ncol = 4) +
    coord_flip() + 
    theme_bw()
  
  out_file <- paste0(out_folder, "/", "barplots.n_top=", n_top, ".", w_col, ".pdf")
  pdf(out_file, width = 8, height = 8)
  plot(plt_up / plt_down)
  dev.off()

}


#~~~~~~~~~~~~~~
#WM-2d centered
###############

# Plots based on WM 2d differences

d_sub_bp <-res[["WM-2d"]]

# Series of bar plots

#ylim <- c(-2.5, 2.5)
#ggplot(d_sub_bp %>% dplyr::filter(name %in% c(sample_w, "Tam-2d", "Fulv-2d")), aes(x=gene, y=value)) + 
#  geom_bar(stat = "identity") +
#  ylim(ylim) +
#  coord_flip() +
#  facet_grid(~name) +
#  theme_bw()

#These bar plots are difficult to grasp
#Try box plots

#ggplot(d_sub_bp %>% dplyr::filter(group != "else"), aes(x=name, y=value)) + 
#  geom_boxplot() + 
#  coord_flip() + 
#  facet_grid(~group) +
#  theme_bw() +
#  ylim(-2.5,2.5) +
#  geom_hline(yintercept=0)

# Pairs (scatters)

d_sub_bp_wide <- d_sub_bp %>% 
  pivot_wider(names_from = name, values_from = value)

#highlight top 10 per up/down
top_10 <- d_sub_bp_wide %>% arrange(`WM-2d`) %>% pull(gene) %>% tail(10) %>% as.character()
bottom_10 <- d_sub_bp_wide %>% arrange(`WM-2d`) %>% pull(gene) %>% head(10) %>% as.character()
d_sub_bp_wide <- d_sub_bp_wide %>% mutate(highlight = gene %in% c(top_10, bottom_10))

out_file <- paste0(out_folder, "/", "data.WM-2d.wide.anno.txt")
write_tsv(x = d_sub_bp_wide, file = out_file)

lim <- c(-3,3)
plts <- list()

#similarities 1st line

plts[["1st_WM_Tam"]] <- customScatter(data = d_sub_bp_wide, 
                                      x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                      y_id = "Tam-2d", y_lab = "Tam-2d vs RM", lim = lim)

plts[["1st_WM_Fulv"]] <- customScatter(data = d_sub_bp_wide, 
                                       x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                       y_id = "Fulv-2d", y_lab = "Fulv-2d vs RM", lim = lim)

#sustained in Resistant

plts[["res_WM_Lted"]] <- customScatter(data = d_sub_bp_wide, 
                                       x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                       y_id = "Lted", y_lab = "Lted vs RM", lim = lim)

plts[["res_Tam_TamR"]] <- customScatter(data = d_sub_bp_wide, 
                                        x_id = "Tam-2d", x_lab = "Tam-2d vs RM", 
                                        y_id = "TamR", y_lab = "TamR vs RM", lim = lim)

plts[["res_Fulv_FulvR"]] <- customScatter(data = d_sub_bp_wide,
                                          x_id = "Fulv-2d", x_lab = "Fulv-2d vs RM", 
                                          y_id = "FulvR", y_lab = "FulvR vs RM", lim = lim)

#1st line (WM) vs 2nd line

plts[["2nd_WM_Tam"]] <- customScatter(data = d_sub_bp_wide, 
                                      x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                      y_id = "Lted_Tam-7d", y_lab = "Lted-Tam-7d vs Lted", lim = lim)

plts[["2nd_WM_Fulv"]] <- customScatter(data = d_sub_bp_wide, 
                                      x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                      y_id = "Lted_Fulv-7d", y_lab = "Lted-Fulv-7d vs Lted", lim = lim)

plts[["2nd_WM_Palbo"]] <- customScatter(data = d_sub_bp_wide, 
                                        x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                        y_id = "Lted_Palbo-7d", y_lab = "Lted-Palbo-7d vs Lted", lim = lim)

plts[["2nd_WM_CDK7i"]] <- customScatter(data = d_sub_bp_wide, 
                                        x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                        y_id = "Lted_CDK7i-7d", y_lab = "Lted-CDK7i-7d vs Lted", lim = lim)

#1st line (Tam) vs 2nd line

plts[["2nd_Tam_WM"]] <- customScatter(data = d_sub_bp_wide, 
                                      x_id = "Tam-2d", x_lab = "Tam-2d vs RM", 
                                      y_id = "TamR_WM-7d", y_lab = "TamR-WM-7d vs TamR", lim = lim)

plts[["2nd_Tam_Fulv"]] <- customScatter(data = d_sub_bp_wide, 
                                      x_id = "Tam-2d", x_lab = "Tam-2d vs RM", 
                                      y_id = "TamR_Fulv-7d", y_lab = "TamR-Fulv-7d vs TamR", lim = lim)

plts[["2nd_Tam_Palbo"]] <- customScatter(data = d_sub_bp_wide, 
                                         x_id = "Tam-2d", x_lab = "Tam-2d vs RM", 
                                         y_id = "TamR_Palbo-7d", y_lab = "TamR-Palbo-7d vs TamR", lim = lim)

plts[["2nd_Tam_CDK7i"]] <- customScatter(data = d_sub_bp_wide, 
                                         x_id = "Tam-2d", x_lab = "Tam-2d vs RM", 
                                         y_id = "TamR_CDK7i-7d", y_lab = "TamR-CDK7i-7d vs TamR", lim = lim)

#Mutant vs wt - 1st line

plts[["mut_RM"]] <- customScatter(data = d_sub_bp_wide, 
                                  x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                  y_id = "D538G_RM", y_lab = "D538G-RM vs RM", lim = lim)

plts[["mut_WM"]] <- customScatter(data = d_sub_bp_wide, 
                                  x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                  y_id = "D538G_WM-2d", y_lab = "D538G-WM vs D538G-RM", lim = lim)

#Neo-adjuvant up/down TF by activity

d_sub_bp_wide_neo <- d_sub_bp_wide %>% 
  select(c("gene", "group", "WM-2d", "Tam-2d", "Fulv-2d", "highlight")) %>%
  inner_join(d_neo, by = "gene")

plts[["neo_WM"]] <- customScatter(data = d_sub_bp_wide_neo, 
                                  x_id = "WM-2d", x_lab = "WM-2d vs RM", 
                                  y_id = "score", y_lab = "Neo-ET Matched", lim = c(-5,9))

#plots
out_file <- paste0(out_folder, "/", "scatter.WM-2d.first_line.pdf")
pdf(out_file, width = 8, height = 4)
w <- grep("1st", names(plts))
for (plt_id in names(plts)[w]) {
  plot(plts[[plt_id]][[1]] + plts[[plt_id]][[2]])
}
dev.off()
out_file <- paste0(out_folder, "/", "scatter.WM-2d.second_line.pdf")
pdf(out_file, width = 8, height = 4)
w <- grep("2nd", names(plts))
for (plt_id in names(plts)[w]) {
  plot(plts[[plt_id]][[1]] + plts[[plt_id]][[2]])
}
dev.off()
out_file <- paste0(out_folder, "/", "scatter.WM-2d.resistant.pdf")
pdf(out_file, width = 8, height = 4)
w <- grep("res", names(plts))
for (plt_id in names(plts)[w]) {
  plot(plts[[plt_id]][[1]] + plts[[plt_id]][[2]])
}
dev.off()
out_file <- paste0(out_folder, "/", "scatter.WM-2d.mutant.pdf")
pdf(out_file, width = 8, height = 4)
w <- grep("mut", names(plts))
for (plt_id in names(plts)[w]) {
  plot(plts[[plt_id]][[1]] + plts[[plt_id]][[2]])
}
dev.off()

#DimPlot(so_subset)
#FeaturePlot(so_subset, features = c("AIP"), min.cutoff = -2, max.cutoff = 2) + 
#  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')
#VlnPlot(so_subset, features = c("AIP"), pt.size = 0)

#AIP -> also called XAP2, could be a negative regulator of ERa signaling
#    -> simply anti-correlated to ESR1?
#AHR -> exactly the opposite pattern


#~~~~~~~~~
#Heat maps
##########

#prepare lists of sample ids
ids <- list()
#first line TFs
ids[["First_line"]] <- metadata_ordered %>% dplyr::filter(class == "first_line") %>% pull(cell_type_2) %>% unique()
#prepare mutant TFs
ids[["ESR1_mutant"]] <- "D538G_RM"
#prepare second line TFs
ids[["Second_line"]] <- metadata_ordered %>% dplyr::filter(class == "second_line" & grepl("7d", cell_type_2)) %>% pull(cell_type_2) %>% unique()
#prepare second line TFs, no CDK7i
w <- grepl("CDK7i", ids[["Second_line"]])
ids[["Second_line_no-CDK7i"]] <- ids[["Second_line"]][!w]
#prepare second line TFs, CDK7i only
ids[["Second_line_only-CDK7i"]] <- ids[["Second_line"]][w]

#add data to res
for (id in names(ids)) {
  any_up <- res_up[names(res_up) %in% ids[[id]]] %>% unlist() %>% as.character() %>% unique()
  any_down <- res_down[names(res_down) %in% ids[[id]]] %>% unlist() %>% as.character() %>% unique()
  res[[id]] <- d_sub_long %>%
    mutate(group = ifelse(gene %in% any_up, "up", ifelse(gene %in% any_down, "down", "else")))
}

#single samples to be considered
ids_single_sample <- c("WM-2d", "Tam-2d", "Fulv-2d")

#consolidate list of ids
ids_to_consider <- c(ids_single_sample, names(ids))

hm_row_ks <- list()
hm_row_ks[["WM-2d"]] <- 4
hm_row_ks[["Tam-2d"]] <- 4
hm_row_ks[["Fulv-2d"]] <- 4
hm_row_ks[["First_line"]] <- 6
hm_row_ks[["ESR1_mutant"]] <- 3
hm_row_ks[["Second_line"]] <- 5
hm_row_ks[["Second_line_no-CDK7i"]] <- 3
hm_row_ks[["Second_line_only-CDK7i"]] <- 3

hm_heights <- list()
hm_heights[["WM-2d"]] <- 11
hm_heights[["Tam-2d"]] <- 7
hm_heights[["Fulv-2d"]] <- 11
hm_heights[["First_line"]] <- 16
hm_heights[["ESR1_mutant"]] <- 6
hm_heights[["Second_line"]] <- 14
hm_heights[["Second_line_no-CDK7i"]] <- 13
hm_heights[["Second_line_only-CDK7i"]] <- 12

d_sub_ct <- c("WM-2d", "Tam-2d", "Fulv-2d",
              "Lted", "Lted_Tam-7d", "Lted_Fulv-7d", "Lted_Palbo-7d", "Lted_CDK7i-7d",
              "TamR", "TamR_WM-7d", "TamR_Fulv-7d", "TamR_Palbo-7d", "TamR_CDK7i-7d",
              "FulvR",
              "D538G_RM", "D538G_WM-2d")

d_sub_ct_groups <- c(rep("first_line", 3), "FulvR", rep("Lted", 5), rep("TamR", 5), rep("Mutant", 2))
d_sub_ct_groups <- factor(d_sub_ct_groups, levels = c("first_line", "Lted", "TamR", "FulvR", "Mutant"))

for (id in ids_to_consider) {
  
  d_sub_bp <- res[[id]]
  
  rows_f <- rownames(d_sub) %in% c(d_sub_bp %>% dplyr::filter(group != "else") %>% pull(gene) %>% as.character())
  cols_order <- d_sub_ct
  d_sub_hm <- d_sub[rows_f,cols_order]
  
  neo_up <- d_neo %>% dplyr::filter(score > 0 & p_value_BH <= 0.05) %>% pull(gene)
  neo_down <- d_neo %>% dplyr::filter(score < 0 & p_value_BH <= 0.05) %>% pull(gene)
  row_anno_neo <- rep("nc", nrow(d_sub_hm))
  row_anno_neo[rownames(d_sub_hm) %in% neo_up] <- "up"
  row_anno_neo[rownames(d_sub_hm) %in% neo_down] <- "down"
  
  ha_cols <- list(NeoAdj_ET = c("up" = "orange", "down" = "darkgreen", "nc" = "lightgrey"))
  
  ha <- rowAnnotation(NeoAdj_ET = row_anno_neo,
                      annotation_name_gp = gpar(fontsize = 8),
                      gap = unit(1, "mm"),
                      col = ha_cols)
  
  d_sub_hm[is.na(d_sub_hm)] <- 0
  
  set.seed(123)
  HM <- Heatmap(d_sub_hm,
                heatmap_legend_param = list(title = "Regulator Activity\n(vs baseline)"),
                col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red")),
                row_names_gp = grid::gpar(fontsize = 8),
                column_title_rot = 90,
                row_title_rot = 0,
                cluster_columns = FALSE,
                column_split= d_sub_ct_groups,
                right_annotation = ha,
                row_km = hm_row_ks[[id]],
                border = TRUE)
  
  out_file <- paste0(out_folder, "/heatmap.", id, ".pdf")
  pdf(out_file, width = 6, height = hm_heights[[id]])
  draw(HM)
  dev.off()
  
  #box plots
  
  plt_bp <- d_sub_bp %>% 
    dplyr::filter(group != "else") %>% 
    mutate(name = factor(name, levels = d_sub_ct)) %>%
    ggplot(aes(x=name, y=abs(value))) + 
    geom_boxplot() + 
    facet_wrap(~group) + 
    coord_flip() + 
    ylim(c(0,2)) + 
    xlab("") +
    ylab("Regulator Activity (vs baseline) [absolute]") +
    theme_bw()
  
  out_file <- paste0(out_folder, "/boxplots.", id, ".pdf")
  pdf(out_file, width = 5, height = 4)
  plot(plt_bp)
  dev.off()
  
}


#~~~~~~~~~~~
#Other plots
############

## mutant ~ first_line, color-coded by neo-adjuvant

plt <- d_sub_wide %>% 
  dplyr::filter(gene %in% tf_w) %>%
  mutate(gene_label = ifelse(abs(`D538G_RM`) >= 1 | abs(first_line) >= 1, gene, NA)) %>% 
  ggplot(aes(x=first_line, y=`D538G_RM`, color=neoadj_class, label=gene_label)) + 
  geom_point() +
  geom_text_repel(max.overlaps = 30) +
  geom_hline(yintercept=0, linetype=2, color = "grey") +
  geom_vline(xintercept=0, linetype=2, color = "grey") +
  scale_color_manual(values=c('darkgreen', 'darkgrey', 'orange')) +
  xlim(c(-3, 3)) +
  ylim(c(-5, 10)) +
  theme_bw()

out_file <- paste0(out_folder, "/", "scatter.first_line.mutant.pdf")
pdf(out_file, width = 7, height = 5.5)
plot(plt)
dev.off()

## hclust of samples

#filter for TFs expressed & active 
d_sub_f <- d_sub[tf_w,]

hc_res <- hclust(dist(t(d_sub_f)))

out_file <- paste0(out_folder, "/", "hclust.samples.pdf")
pdf(out_file, width = 8, height = 6)
plot(hc_res)
dev.off()

## heat map centered on neo-adjuvant TFs alone (flagged by machine learning)

neo_up <- d_neo %>% dplyr::filter(score > 0 & p_value_BH <= 0.05) %>% pull(gene)
neo_down <- d_neo %>% dplyr::filter(score < 0 & p_value_BH <= 0.05) %>% pull(gene)

neo_sel_lasso <- d_neo %>% dplyr::filter(lasso_est != 0) %>% pull(gene)
neo_top_rf <- d_neo %>% arrange(-rf_imp) %>% head(20) %>% pull(gene)

rows_f <- rownames(d_sub_f) %in% c(neo_up, neo_down)
cols_order <- d_sub_ct
d_sub_hm <- d_sub_f[rows_f,cols_order]

row_anno_neo <- rep("nc", nrow(d_sub_hm))
row_anno_neo[rownames(d_sub_hm) %in% neo_up] <- "up"
row_anno_neo[rownames(d_sub_hm) %in% neo_down] <- "down"

row_anno_lasso <- rep("no", nrow(d_sub_hm))
row_anno_lasso[rownames(d_sub_hm) %in% neo_sel_lasso] <- "yes"

row_anno_rf <- rep("no", nrow(d_sub_hm))
row_anno_rf[rownames(d_sub_hm) %in% neo_top_rf] <- "yes"

ha_cols <- list(NeoAdj_ET = c("up" = "orange", "down" = "darkgreen", "nc" = "lightgrey"),
                NeoAdj_Lasso = c("yes" = "blue", "no" = "lightgrey"),
                NeoAdj_RF = c("yes" = "navy", "no" = "lightgrey"))

ha <- rowAnnotation(NeoAdj_ET = row_anno_neo,
                    NeoAdj_Lasso = row_anno_lasso,
                    NeoAdj_RF = row_anno_rf,
                    annotation_name_gp = gpar(fontsize = 8),
                    gap = unit(1, "mm"),
                    col = ha_cols)

d_sub_hm[is.na(d_sub_hm)] <- 0

set.seed(123)
HM <- Heatmap(d_sub_hm,
              heatmap_legend_param = list(title = "Regulator Activity\n(vs baseline)"),
              col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red")),
              row_names_gp = grid::gpar(fontsize = 8),
              column_title_rot = 90,
              row_title_rot = 0,
              cluster_columns = FALSE,
              column_split= d_sub_ct_groups,
              right_annotation = ha,
              border = TRUE)

out_file <- paste0(out_folder, "/heatmap.NeoAdj.pdf")
pdf(out_file, width = 5.5, height = 10)
draw(HM)
dev.off()

