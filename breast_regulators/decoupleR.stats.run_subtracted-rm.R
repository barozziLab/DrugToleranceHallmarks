
out_folder <- "decoupleR.stats.run_subtracted-rm/"
cmd <- paste0("mkdir ", out_folder)
system(cmd)


#~~~
#Any
####

# Analyses based on differences in any of the conditions

cols_w <- d_sub_rm_long$name %>% unique() %>% as.character()

res <- list()
any_up <- c()
any_down <- c()

for (col_w in cols_w) {
  
  res[[col_w]] <- data_filter(d_sub_long = d_sub_rm_long, 
                              tf_expr_frac = tf_expr_frac, 
                              sample_w = col_w, 
                              filt_tf_expr = filt_tf_expr,
                              filt_tf_act = filt_tf_act)
  
  reg_up <- res[[col_w]] %>% dplyr::filter(group == "up") %>% pull(gene) %>% unique() %>% as.character()
  any_up <- c(any_up, reg_up)
  
  reg_down <- res[[col_w]] %>% dplyr::filter(group == "down") %>% pull(gene) %>% unique() %>% as.character()
  any_down <- c(any_down, reg_down)
  
}

out_file <- paste0(out_folder, "/", "data.all.res.rds")
saveRDS(object = res, file = out_file)

any_up <- unique(sort(any_up))
any_down <- unique(sort(any_down))

#Add the grouping

d_sub_bp <- d_sub_rm_long
d_sub_bp <- d_sub_bp %>%
  mutate(group = ifelse(gene %in% any_up, "up", ifelse(gene %in% any_down, "down", "else")))


# Heat map

d_sub_ct <- c("WM-2d", "Tam-2d", "Fulv-2d",
              "Lted", "Lted_Tam-7d", "Lted_Fulv-7d", "Lted_Palbo-7d", "Lted_CDK7i-7d",
              "TamR", "TamR_WM-7d", "TamR_Fulv-7d", "TamR_Palbo-7d", "TamR_CDK7i-7d",
              "FulvR",
              "D538G_RM", "D538G_WM-2d")

d_sub_ct_groups <- c(rep("first_line", 3), "FulvR", rep("Lted", 5), rep("TamR", 5), rep("Mutant", 2))
d_sub_ct_groups <- factor(d_sub_ct_groups, levels = c("first_line", "Lted", "TamR", "FulvR", "Mutant"))

rows_f <- rownames(d_sub_rm) %in% c(d_sub_bp %>% dplyr::filter(group != "else") %>% pull(gene) %>% as.character())
cols_order <- d_sub_ct
d_sub_hm <- d_sub_rm[rows_f,cols_order]

neo_up <- d_neo %>% dplyr::filter(score > 0 & p_value_BH <= 0.05) %>% pull(gene)
neo_down <- d_neo %>% dplyr::filter(score < 0 & p_value_BH <= 0.05) %>% pull(gene)
row_anno_neo <- rep("nc", nrow(d_sub_hm))
row_anno_neo[rownames(d_sub_hm) %in% neo_up] <- "up"
row_anno_neo[rownames(d_sub_hm) %in% neo_down] <- "down"

ha_cols <- list(NeoAdj_ET = c("up" = "orange", "down" = "darkgreen", "nc" = "lightgrey"))

ha <- rowAnnotation(NeoAdj_ET = row_anno_neo,
                    annotation_name_gp= gpar(fontsize = 8),
                    gap = unit(1, "mm"),
                    col = ha_cols)

set.seed(123)
HM <- Heatmap(d_sub_hm,
              heatmap_legend_param = list(title = "Regulator Activity\n(vs RM)"),
              col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red")),
              row_names_gp = grid::gpar(fontsize = 0),
              column_title_rot = 90,
              row_title_rot = 0,
              cluster_columns = FALSE,
              column_split= d_sub_ct_groups,
              right_annotation = ha,
              border = TRUE)

out_file <- paste0(out_folder, "/", "heatmap.any.compact.pdf")
pdf(out_file, width = 5, height = 8)
draw(HM)
dev.off()

set.seed(123)
HM <- Heatmap(d_sub_hm,
              heatmap_legend_param = list(title = "Regulator Activity\n(vs RM)"),
              col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red")),
              row_names_gp = grid::gpar(fontsize = 8),
              column_title_rot = 90,
              row_title_rot = 0,
              cluster_columns = FALSE,
              column_split= d_sub_ct_groups,
              right_annotation = ha,
              border = TRUE)

out_file <- paste0(out_folder, "/", "heatmap.any.verbose.pdf")
pdf(out_file, width = 5.5, height = 19)
draw(HM)
dev.off()


#~~~~~~~~~~~
#Other plots
############

## hclust of samples

#filter for TFs expressed & active 
d_sub_rm_f <- d_sub_rm[tf_w,]

hc_res <- hclust(dist(t(d_sub_rm_f)))

out_file <- paste0(out_folder, "/", "hclust.samples.pdf")
pdf(out_file, width = 8, height = 5)
plot(hc_res)
dev.off()

## heat map centered on neo-adjuvant TFs alone (flagged by machine learning)

neo_up <- d_neo %>% dplyr::filter(score > 0 & p_value_BH <= 0.05) %>% pull(gene)
neo_down <- d_neo %>% dplyr::filter(score < 0 & p_value_BH <= 0.05) %>% pull(gene)

neo_sel_lasso <- d_neo %>% dplyr::filter(lasso_est != 0) %>% pull(gene)
neo_top_rf <- d_neo %>% arrange(-rf_imp) %>% head(20) %>% pull(gene)

rows_f <- rownames(d_sub_rm_f) %in% c(neo_up, neo_down)
cols_order <- d_sub_ct
d_sub_hm <- d_sub_rm_f[rows_f,cols_order]

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
pdf(out_file, width = 5.5, height = 11)
draw(HM)
dev.off()

