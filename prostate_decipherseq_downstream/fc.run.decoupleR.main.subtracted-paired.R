
out_folder <- "fc_results/decoupleR.run_subtracted-paired/"
cmd <- paste0("mkdir ", out_folder)
system(cmd)


#~~~~~~~~~~~~~~~~
#Data preparation
#################

# Overall, annotated dataset (TF by sample)

d_sub_wide <- d_sub %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble()

#add average first and second line
d_sub_wide <- d_sub_wide %>% 
  rowwise() %>% 
  mutate(first_line = `LNCaP_WM-7d`) %>%
  mutate(second_line = mean(c(`LNCaP95_Enz-7d`, 
                              `LNCaP95_WM-7d`)))

# Analyses based on differences in any of the conditions

cols_w <- d_sub_long$name %>% unique() %>% as.character()

res <- list()

res_up <- c()
res_down <- c()

for (col_w in cols_w) {
  
  res[[col_w]] <- data_filter(d_sub_long = d_sub_long, 
                              tf_expr_frac = tf_expr_frac$`cell-type`, 
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
ids[["First_line"]] <- metadata_ordered %>% dplyr::filter(class == "first_line") %>% pull(description) %>% unique()
#second line
ids[["Second_line"]] <- metadata_ordered %>% dplyr::filter(class == "second_line") %>% pull(description) %>% unique()
#resistant
ids[["Resistant"]] <- metadata_ordered %>% dplyr::filter(class == "resistant") %>% pull(description) %>% unique()

ids_all <- unlist(ids)  %>% as.character() %>% sort() %>% unique()
ids_all <- ids_all[ids_all != "LNCaP_WM-2d"]

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

#Sorted on First_line, etc

for (w_col in c("First_line", "Second_line", "Resistant")) {
  
  goi <- d_goi_summary_w_verb %>% arrange(-get(w_col)) %>% pull(gene) %>% head(n_top) %>% rev()
  
  plt_up <- d_goi_summary %>% 
    dplyr::filter(gene %in% goi) %>% 
    mutate(gene = factor(gene, levels = goi)) %>%
    ggplot(aes(x=gene, y=mean_value)) +
    geom_bar(stat="identity") + 
    ylim(-2.5, 4) + 
    xlab("Top Regulator (up)") +
    ylab("Mean Activity (vs baseline)") +
    facet_wrap(~sample_cls, ncol = 4) +
    coord_flip() + 
    theme_bw()
  
  goi <- d_goi_summary_w_verb %>% arrange(get(w_col)) %>% pull(gene) %>% head(n_top)
  
  plt_down <- d_goi_summary %>% 
    dplyr::filter(gene %in% goi) %>% 
    mutate(gene = factor(gene, levels = goi)) %>%
    ggplot(aes(x=gene, y=mean_value)) +
    geom_bar(stat="identity") + 
    ylim(-2.5, 4) + 
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
#WM-7d centered
###############

# Plots based on WM 7d differences

d_sub_bp <-res[["LNCaP_WM-7d"]]

# Series of bar plots

#ylim <- c(-2.5, 2.5)
#ggplot(d_sub_bp %>% dplyr::filter(name %in% c(sample_w, "MCF7_Tam-2d", "MCF7_Fulv-2d")), aes(x=gene, y=value)) + 
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
top_10 <- d_sub_bp_wide %>% arrange(`LNCaP_WM-7d`) %>% pull(gene) %>% tail(10) %>% as.character()
bottom_10 <- d_sub_bp_wide %>% arrange(`LNCaP_WM-7d`) %>% pull(gene) %>% head(10) %>% as.character()
d_sub_bp_wide <- d_sub_bp_wide %>% mutate(highlight = gene %in% c(top_10, bottom_10))

out_file <- paste0(out_folder, "/", "data.LNCaP_WM-7d.wide.anno.txt")
write_tsv(x = d_sub_bp_wide, file = out_file)

lim <- c(-3,3)
plts <- list()

#sustained in Resistant

plts[["res_WM_LNCaP95"]] <- customScatter(data = d_sub_bp_wide, 
                                          x_id = "LNCaP_WM-7d", x_lab = "WM-7d vs RM", 
                                          y_id = "LNCaP95", y_lab = "LNCaP95 vs RM", lim = lim)

#1st line (WM) vs 2nd line

plts[["2nd_WM_1st_WM"]] <- customScatter(data = d_sub_bp_wide, 
                                         x_id = "LNCaP_WM-7d", x_lab = "WM-7d vs RM", 
                                         y_id = "LNCaP95_WM-7d", y_lab = "LNCaP95-WM-7d vs LNCaP95", lim = lim)

plts[["2nd_WM_1st_Enz"]] <- customScatter(data = d_sub_bp_wide, 
                                          x_id = "LNCaP_WM-7d", x_lab = "WM-7d vs RM", 
                                          y_id = "LNCaP95_Enz-7d", y_lab = "LNCaP95-Enz-7d vs LNCaP95", lim = lim)

#similarities 2nd line

plts[["2nd_WM_Enz"]] <- customScatter(data = d_sub_bp_wide, 
                                      x_id = "LNCaP95_WM-7d", x_lab = "LNCaP95-WM-7d vs LNCaP95", 
                                      y_id = "LNCaP95_Enz-7d", y_lab = "LNCaP95-Enz-7d vs LNCaP95", lim = lim)

#plots
out_file <- paste0(out_folder, "/", "scatterplots.pdf")
pdf(out_file, width = 8, height = 4)
plot(plts[["res_WM_LNCaP95"]][[1]] + plts[["res_WM_LNCaP95"]][[2]])
plot(plts[["2nd_WM_1st_WM"]][[1]] + plts[["2nd_WM_1st_WM"]][[2]])
plot(plts[["2nd_WM_1st_Enz"]][[1]] + plts[["2nd_WM_1st_Enz"]][[2]])
plot(plts[["2nd_WM_Enz"]][[1]] + plts[["2nd_WM_Enz"]][[2]])
dev.off()

#SREBF1 and 2, up in first line, sustained in resistant, but down in second line

#DimPlot(so)
#FeaturePlot(so, features = c("SREBF1"), min.cutoff = -2, max.cutoff = 2) + 
#  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')
#VlnPlot(so, features = c("SREBF1"), pt.size = 0)

