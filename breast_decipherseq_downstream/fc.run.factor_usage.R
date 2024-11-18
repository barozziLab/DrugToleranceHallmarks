
col_fun <- colorRamp2(breaks = seq(0,10,1), colors = rev(brewer.pal(11,"Spectral")))
col_fun_scaled <- colorRamp2(c(-3, 0, 3), colors = c("navy", "white", "red"))
col_fun_scaled_short <- colorRamp2(c(-1.5, 0, 1.5), colors = c("navy", "white", "red"))


#~~~~~~~~~~~~~~~~
#Prepare Metadata
#################

metadata_ordered <- read.xls("metadata_filtered_ordered.xlsx") %>% 
  as_tibble() %>%
  mutate(sample = paste0(description, "_batch-", batch)) %>%
  mutate(class = ifelse(ESR1_status == "mut", "ESR1_mutant", class)) %>%
  mutate(class = factor(class, levels = c("treatment-naive", "ESR1_mutant", "first_line", "resistant", "second_line")))

metadata_unique <- metadata_ordered %>% 
  select(description, treatment, ESR1_status, class) %>%
  unique() %>%
  mutate(class_acute = ifelse(class %in% c("first_line", "second_line"), TRUE, FALSE)) %>%
  mutate(class_resistant = ifelse(class == "resistant", TRUE, FALSE)) %>%
  mutate(class_resistant_lted = ifelse(grepl("Lted", description), TRUE, FALSE)) %>%
  mutate(class_resistant_tamR = ifelse(grepl("TamR", description), TRUE, FALSE)) %>%
  mutate(class_resistant_fulvR = ifelse(grepl("FulvR", description), TRUE, FALSE)) %>%
  mutate(class_treatment_naive = ifelse(class == "treatment-naive", TRUE, FALSE)) %>%
  mutate(class_first_line = ifelse(class == "first_line", TRUE, FALSE)) %>%
  mutate(class_second_line = ifelse(class == "second_line", TRUE, FALSE)) %>%
  mutate(class_esr1_mutant = ifelse(ESR1_status == "mut", TRUE, FALSE)) %>%
  select(-ESR1_status)

desc <- metadata_unique$description %>% unique()

for (desc_i in desc) {
  field_name <- paste0("condition_", desc_i)
  metadata_unique <- metadata_unique %>% 
    mutate(!!field_name := ifelse(description == desc_i, TRUE, FALSE))
}


#~~~~~~~~~~~~~~~~~~
#Summary Statistics
###################

# stats
prog_stats <- H.data.long %>%
  group_by(sample, gene_program) %>%
  summarize(mean = mean(value)) %>% 
  mutate(gene_program = as.numeric(str_remove(gene_program, "R33_Program"))) %>%
  ungroup()

# mean program usage
mean_prog_usage <- prog_stats %>%
  pivot_wider(names_from = gene_program, values_from = mean) %>%
  column_to_rownames("sample")

# numerically order column names
col_ord <- as.character(sort(as.numeric(colnames(mean_prog_usage))))
mean_prog_usage <- mean_prog_usage[, col_ord] 

# sort samples according to metadata
mean_prog_usage_ordered <- mean_prog_usage[metadata_ordered$sample, ]

hm_mean_linear <- 
  Heatmap(mean_prog_usage_ordered, 
          col = col_fun, 
          cluster_rows = FALSE, 
          name = "mean_program_usage", 
          row_split = metadata_ordered$class,
          row_title_rot = 0, 
          row_gap = unit(3, "mm"),
          heatmap_legend_param = list(title = "mean prog. usage"),
          row_names_max_width = unit(10, "cm"),
          border = TRUE)

cols_order <- column_order(draw(hm_mean_linear))
programs_order <- colnames(mean_prog_usage_ordered)[cols_order]

# scale columns
mean_prog_usage_ordered_scaled <- apply(t(mean_prog_usage_ordered), 1, scale)
rownames(mean_prog_usage_ordered_scaled) <- rownames(mean_prog_usage_ordered)

# note: ordered according to the hierarchical clustering of the linear matrix
hm_mean_scaled <- 
  Heatmap(mean_prog_usage_ordered_scaled[,cols_order], 
          col = col_fun_scaled, 
          cluster_rows = FALSE, 
          cluster_columns = FALSE, 
          name = "mean_program_usage", 
          row_split = metadata_ordered$class,
          row_title_rot = 0, 
          row_gap = unit(3, "mm"),
          heatmap_legend_param = list(title = "mean prog. usage"),
          row_names_max_width = unit(10, "cm"),
          border = TRUE)

# subtract the relevant baseline

mean_prog_usage_ordered_sub <- mean_prog_usage_ordered
#
#first line - rm [baseline = rm average] 
w <- metadata_ordered %>% dplyr::filter(class == "treatment-naive") %>% pull(sample)
baseline <- apply(mean_prog_usage_ordered[w,], 2, mean)
w <- metadata_ordered %>% dplyr::filter(class == "treatment-naive" | class == "first_line") %>% pull(sample)
mean_prog_usage_ordered_sub[w,] <- t(t(mean_prog_usage_ordered[w,]) - baseline)
#
#mutant wm - mutant rm
mean_prog_usage_ordered_sub["MCF7-D538G_RM_batch-2020B",] <- mean_prog_usage_ordered["MCF7-D538G_RM_batch-2020B",] - mean_prog_usage_ordered["MCF7-D538G_RM_batch-2020B",]
mean_prog_usage_ordered_sub["MCF7-D538G_WM-2d_batch-2020B",] <- mean_prog_usage_ordered["MCF7-D538G_WM-2d_batch-2020B",] - mean_prog_usage_ordered["MCF7-D538G_RM_batch-2020B",]
#
#resistant - rm [baseline = rm average]
w <- metadata_ordered %>% dplyr::filter(class == "treatment-naive") %>% pull(sample)
baseline <- apply(mean_prog_usage_ordered[w,], 2, mean)
w <- metadata_ordered %>% dplyr::filter(class == "treatment-naive" | class == "resistant") %>% pull(sample)
mean_prog_usage_ordered_sub[w,] <- t(t(mean_prog_usage_ordered[w,]) - baseline)
#
#second line - resistant [baseline = resistant average, separate for each drug]
w <- metadata_ordered %>% dplyr::filter(grepl("TamR_batch", sample)) %>% pull(sample)
baseline <- apply(mean_prog_usage_ordered[w,], 2, mean)
w <- metadata_ordered %>% dplyr::filter(grepl("TamR", description)) %>% pull(sample)
mean_prog_usage_ordered_sub[w,] <- t(t(mean_prog_usage_ordered[w,]) - baseline)
#
w <- metadata_ordered %>% dplyr::filter(grepl("Lted_batch", sample)) %>% pull(sample)
baseline <- apply(mean_prog_usage_ordered[w,], 2, mean)
w <- metadata_ordered %>% dplyr::filter(grepl("Lted", description)) %>% pull(sample)
mean_prog_usage_ordered_sub[w,] <- t(t(mean_prog_usage_ordered[w,]) - baseline)

# note: ordered according to the hierarchical clustering of the linear matrix
hm_mean_sub <- 
  Heatmap(mean_prog_usage_ordered_sub[,cols_order], 
          col = col_fun_scaled, 
          cluster_rows = FALSE, 
          cluster_columns = FALSE, 
          name = "mean_program_usage", 
          row_split = metadata_ordered$class,
          row_title_rot = 0, 
          row_gap = unit(3, "mm"),
          heatmap_legend_param = list(title = "mean prog. usage"),
          row_names_max_width = unit(10, "cm"),
          border = TRUE)

#plt_mean_ridge <- ggplot(prog_stats %>% mutate(gene_program = as.factor(gene_program)), aes(x = mean, y = gene_program)) + 
#  geom_density_ridges() +
#  scale_x_log10() +
#  xlab("Program Activity (log10)") +
#  ylab("Gene Program") +
#  theme_bw()

gene_program_levels <- colnames(mean_prog_usage_ordered)[cols_order]
plt_mean_bplot <- ggplot(prog_stats %>% mutate(gene_program = factor(gene_program, levels = gene_program_levels)), aes(x = gene_program, y = mean)) + 
  geom_boxplot() +
  ylab("Program Activity") +
  xlab("Gene Program") +
  ylim(0,10) +
  theme_bw()
plt_mean_bplot_log10 <- ggplot(prog_stats %>% mutate(gene_program = factor(gene_program, levels = gene_program_levels)), aes(x = gene_program, y = mean)) + 
  geom_boxplot() +
  scale_y_log10() +
  ylab("Program Activity (log10)") +
  xlab("Gene Program") +
  theme_bw()

out_path <- paste0(out_folder, "/H.mean.heatmap.linear.pdf")  
pdf(out_path, width = 12, height = 8)
draw(hm_mean_linear)
dev.off()

out_path <- paste0(out_folder, "/H.mean.heatmap.scaled.pdf")  
pdf(out_path, width = 12, height = 8)
draw(hm_mean_scaled)
dev.off()

out_path <- paste0(out_folder, "/H.mean.heatmap.subtracted.pdf")  
pdf(out_path, width = 12, height = 8)
draw(hm_mean_sub)
dev.off()

out_path <- paste0(out_folder, "/H.mean.boxplots.linear.pdf")  
pdf(out_path, width = 7.5, height = 3)
plot(plt_mean_bplot)
dev.off()

out_path <- paste0(out_folder, "/H.mean.boxplots.log10.pdf")  
pdf(out_path, width = 7.5, height = 3)
plot(plt_mean_bplot_log10)
dev.off()

H.stats <- H.data.long %>%
  select(batch, kit, description, gene_program, value) %>%
  left_join(metadata_unique, by = "description") %>%
  dplyr::rename(program = gene_program, program_value = value, condition = description) %>%
  mutate(program = as.numeric(str_remove(program, "R33_Program")))


#~~~~~~~~~~~~~~~~~~~~
#Program Associations
#####################

## Mutual information + AUC associations

coi <- colnames(H.stats)
coi <- coi[!coi %in% c("program", "program_value")]

coi_auc <- coi[!coi %in% c("batch", "kit", "condition", "treatment", "class")]

stats_res <- c()

program_is <- H.stats$program %>% unique()

for (program_i in program_is) {
  
  input_data <- H.stats %>%
    filter(program == program_i)
  
  for (coi_i in coi) {
    
    #dat <- discretize(input_data[,c("program_value", coi_i)])
    #mi <- mutinformation(dat[,"program_value"], dat[,coi_i], method="emp")
    #stats_res <- rbind(stats_res, c(program_i, coi_i, "MI", mi))
    
    if (coi_i %in% coi_auc) {
      
      res_roc <- roc(input_data[,coi_i,drop=TRUE], input_data[,"program_value",drop=TRUE])
      res_auc <- auc(res_roc)
      stats_res <- rbind(stats_res, c(program_i, coi_i, "AUC", res_auc))
      
    }
    
  }
  
}

stats_res_tib <- tibble(program = stats_res[,1],
                        feature = stats_res[,2],
                        measurement = stats_res[,3],
                        value = as.numeric(stats_res[,4]))

out_path <- paste0(out_folder, "/H.stats.txt") 
write_tsv(x = stats_res_tib, file = out_path)

## Visualize results for AUC (class)

data <- stats_res_tib %>% 
  dplyr::filter(measurement == "AUC" & grepl("class_", feature)) %>%
  dplyr::filter(program %in% programs_order)

lvls_feature <- c("class_treatment_naive", "class_acute", "class_first_line", "class_second_line",
                  "class_resistant", "class_resistant_lted", "class_resistant_tamR", "class_resistant_fulvR",
                  "class_esr1_mutant")

data <- data %>% 
  mutate(value = ifelse(value <= 0.5, 0.5, value)) %>%
  mutate(value = ifelse(value >= 0.9, 0.9, value)) %>%
  mutate(feature = factor(feature, levels = rev(lvls_feature))) %>%
  mutate(program = factor(program, levels = programs_order))

#data_wide <- data %>% 
#  dplyr::select(program, feature, value) %>% 
#  pivot_wider(names_from = feature, values_from = value) %>%
#  column_to_rownames(var = "program")
#data_wide[is.na(data_wide)] <- 0

#hc_rows <- hclust(as.dist(1 - cor(t(data_wide), method = "spearman")))
#hc_cols <- hclust(as.dist(1 - cor(data_wide, method = "spearman")))

plt <- ggplot(data, aes(x = program, y = feature)) + 
  geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "white", high = "red") + 
  scale_size(range = c(.1, 5), name="AUC") +
  theme(legend.key.size = unit(0.4, 'cm'))

out_path <- paste0(out_folder, "/H.stats.class.AUC.bubbleplot.pdf")
pdf(out_path, width = 7.5, height = 2)
plot(plt)
dev.off()

## Visualize results for AUC (condition)

data <- stats_res_tib %>% 
  dplyr::filter(measurement == "AUC" & grepl("condition_", feature)) %>%
  dplyr::filter(program %in% programs_order)

lvls_feature <- data$feature %>% unique()

data <- data %>% 
  mutate(value = ifelse(value <= 0.5, 0.5, value)) %>%
  mutate(value = ifelse(value >= 0.9, 0.9, value)) %>%
  mutate(feature = factor(feature, levels = rev(lvls_feature))) %>%
  mutate(program = factor(program, levels = programs_order))

#data_wide <- data %>% 
#  dplyr::select(program, feature, value) %>% 
#  pivot_wider(names_from = feature, values_from = value) %>%
#  column_to_rownames(var = "program")
#data_wide[is.na(data_wide)] <- 0

#hc_rows <- hclust(as.dist(1 - cor(t(data_wide), method = "spearman")))
#hc_cols <- hclust(as.dist(1 - cor(data_wide, method = "spearman")))

plt <- ggplot(data, aes(x = program, y = feature)) + 
  geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "white", high = "red") + 
  scale_size(range = c(.1, 5), name="AUC") +
  theme(legend.key.size = unit(0.4, 'cm'))

out_path <- paste0(out_folder, "/H.stats.condition.AUC.bubbleplot.pdf")
pdf(out_path, width = 7.5, height = 5)
plot(plt)
dev.off()

plt <- ggplot(data %>% dplyr::filter(!grepl("Lted_.*-2d|TamR_.*-2d", feature)), aes(x = program, y = feature)) + 
  geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "white", high = "red") + 
  scale_size(range = c(.1, 5), name="AUC") +
  theme(legend.key.size = unit(0.4, 'cm'))

out_path <- paste0(out_folder, "/H.stats.condition.AUC.bubbleplot.no_2L-2d.pdf")
pdf(out_path, width = 7.5, height = 3)
plot(plt)
dev.off()

## Visualize results for MI (other features)

#data <- stats_res_tib %>% 
#  dplyr::filter(measurement == "MI" & !grepl("condition_|class_", feature)) %>%
#  dplyr::filter(program %in% programs_order)

#lvls_feature <- data$feature %>% unique()

#data <- data %>% 
#  mutate(feature = factor(feature, levels = rev(lvls_feature))) %>%
#  mutate(program = factor(program, levels = programs_order))

#data_wide <- data %>% 
#  dplyr::select(program, feature, value) %>% 
#  pivot_wider(names_from = feature, values_from = value) %>%
#  column_to_rownames(var = "program")
#data_wide[is.na(data_wide)] <- 0

#hc_rows <- hclust(as.dist(1 - cor(t(data_wide), method = "spearman")))
#hc_cols <- hclust(as.dist(1 - cor(data_wide, method = "spearman")))

#plt <- ggplot(data, aes(x = program, y = feature)) + 
#  geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21) + 
#  theme_bw() + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#  scale_fill_gradient(low = "white", high = "red") + 
#  scale_size(range = c(.1, 5), name="MI") +
#  theme(legend.key.size = unit(0.4, 'cm'))

#out_path <- paste0(out_folder, "/H.stats.global.MI.bubbleplot.pdf")
#pdf(out_path, width = 7.5, height = 2)
#plot(plt)
#dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~
#First-line vs Second-line
##########################

## Subtracted values by group

mean_prog_usage_sub_long <- mean_prog_usage_ordered_sub %>% 
  rownames_to_column(var = "sample") %>% 
  as_tibble() %>%
  pivot_longer(cols = -c("sample"), names_to = "program")

anno <- H.data.long %>% 
  select(sample, description) %>% 
  unique()

mean_prog_usage_sub_long <- mean_prog_usage_sub_long %>%
  left_join(anno %>% unique(), by = "sample") %>%
  left_join(metadata_unique, by = "description")

mean_prog_usage_sub_long$class_global <- "TN"
mean_prog_usage_sub_long$class_global[mean_prog_usage_sub_long$class_first_line] <- "First-Line"
mean_prog_usage_sub_long$class_global[mean_prog_usage_sub_long$class_second_line] <- "Second-Line"
mean_prog_usage_sub_long$class_global[mean_prog_usage_sub_long$class_resistant] <- "Resistant"
mean_prog_usage_sub_long$class_global[mean_prog_usage_sub_long$class_esr1_mutant] <- "ER-mutant"

mean_prog_usage_sub_long$class_global <- factor(mean_prog_usage_sub_long$class_global,
                                                levels = c("TN", "First-Line", "Resistant", "Second-Line", "ER-mutant"))

#sorted programs based on first line median values
programs_sorted_FL <- mean_prog_usage_sub_long %>% 
  group_by(program, class_global) %>% 
  summarise(program_median = median(value)) %>% 
  filter(class_global == "First-Line") %>% 
  arrange(-program_median) %>% 
  pull(program)

mean_prog_usage_sub_long$program <- factor(mean_prog_usage_sub_long$program,
                                           levels = programs_sorted_FL)

#no 2d second-line
plt_mean_bplot_sub <- ggplot(mean_prog_usage_sub_long %>% filter(!(class_second_line & grepl("2d", sample))), aes(x=program, y=value, fill=class_global)) +
  geom_boxplot() +
  ylim(-10,10) +
  theme_bw()

out_path <- paste0(out_folder, "/H.mean.boxplots.subtracted.FL-sorted.pdf")
pdf(out_path, width = 10, height = 4)
plot(plt_mean_bplot_sub)
dev.off()

#median of the box plots, as heat map
mean_prog_usage_sub_mat <- mean_prog_usage_sub_long %>% 
  group_by(class_global, program) %>%
  summarise(value = median(value)) %>%
  pivot_wider(names_from = program, values_from = value) %>%
  column_to_rownames(var = "class_global")

#hm_mean_bplot_sub <- 
#  Heatmap(mean_prog_usage_sub_mat[rownames(mean_prog_usage_sub_mat) != "TN",], 
#          col = col_fun_scaled_short, 
#          cluster_rows = FALSE, 
#          cluster_columns = FALSE, 
#          name = "mean_program_usage", 
#          row_title_rot = 0, 
#          row_gap = unit(3, "mm"),
#          heatmap_legend_param = list(title = "mean prog. usage"),
#          row_names_max_width = unit(10, "cm"),
#          border = TRUE)

#re-plot subtracted heatmap, based on this order, and excluding the 2d second-line
w_samples <- mean_prog_usage_sub_long %>% filter(!(class_second_line & grepl("2d", sample))) %>% pull(sample) %>% unique()
w_rows <- rownames(mean_prog_usage_ordered_sub) %in% w_samples
hm_mean_sub_reord <- 
  Heatmap(mean_prog_usage_ordered_sub[w_rows,programs_sorted_FL], 
          col = col_fun_scaled_short, 
          cluster_rows = FALSE, 
          cluster_columns = FALSE, 
          name = "mean_program_usage", 
          row_split = metadata_ordered$class[w_rows],
          row_title_rot = 0, 
          row_gap = unit(3, "mm"),
          heatmap_legend_param = list(title = "mean prog. usage"),
          row_names_max_width = unit(10, "cm"),
          border = TRUE)

out_path <- paste0(out_folder, "/H.mean.heatmap.subtracted.FL-sorted.pdf")  
pdf(out_path, width = 12, height = 6)
draw(hm_mean_sub_reord)
dev.off()

#AUC plots (class)

data <- stats_res_tib %>% 
  dplyr::filter(measurement == "AUC" & grepl("class_", feature)) %>%
  dplyr::filter(program %in% programs_order)

lvls_feature <- c("class_treatment_naive", "class_acute", "class_first_line", "class_second_line",
                  "class_resistant", "class_resistant_lted", "class_resistant_tamR", "class_resistant_fulvR",
                  "class_esr1_mutant")

data <- data %>% 
  mutate(value = ifelse(value <= 0.5, 0.5, value)) %>%
  mutate(value = ifelse(value >= 0.9, 0.9, value)) %>%
  mutate(feature = factor(feature, levels = rev(lvls_feature))) %>%
  mutate(program = factor(program, levels = programs_sorted_FL))

plt <- ggplot(data, aes(x = program, y = feature)) + 
  geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "white", high = "red") + 
  scale_size(range = c(.1, 5), name="AUC") +
  theme(legend.key.size = unit(0.4, 'cm'))

out_path <- paste0(out_folder, "/H.stats.class.AUC.bubbleplot.FL-sorted.pdf")
pdf(out_path, width = 7.5, height = 2)
plot(plt)
dev.off()

#AUC plots (condition)

data <- stats_res_tib %>% 
  dplyr::filter(measurement == "AUC" & grepl("condition_", feature)) %>%
  dplyr::filter(program %in% programs_order)

lvls_feature <- data$feature %>% unique()

data <- data %>% 
  mutate(value = ifelse(value <= 0.5, 0.5, value)) %>%
  mutate(value = ifelse(value >= 0.9, 0.9, value)) %>%
  mutate(feature = factor(feature, levels = rev(lvls_feature))) %>%
  mutate(program = factor(program, levels = programs_sorted_FL))

plt <- ggplot(data, aes(x = program, y = feature)) + 
  geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "white", high = "red") + 
  scale_size(range = c(.1, 5), name="AUC") +
  theme(legend.key.size = unit(0.4, 'cm'))

out_path <- paste0(out_folder, "/H.stats.condition.AUC.bubbleplot.FL-sorted.pdf")
pdf(out_path, width = 7.5, height = 5)
plot(plt)
dev.off()

plt <- ggplot(data %>% dplyr::filter(!grepl("Lted_.*-2d|TamR_.*-2d", feature)), aes(x = program, y = feature)) + 
  geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "white", high = "red") + 
  scale_size(range = c(.1, 5), name="AUC") +
  theme(legend.key.size = unit(0.4, 'cm'))

out_path <- paste0(out_folder, "/H.stats.condition.AUC.bubbleplot.FL-sorted.no_2L-2d.pdf")
pdf(out_path, width = 7.5, height = 3)
plot(plt)
dev.off()

## Consider only those programs with high AUC for first or second line conditions

#stats_res_tib_poi <- stats_res_tib %>% 
#  dplyr::filter(measurement == "AUC" & grepl("^condition_", feature)) %>% 
#  dplyr::filter(value >= 0.8) %>% rowwise() %>% 
#  mutate(feature = gsub("condition_", "", feature)) %>% 
#  ungroup()
#p_all <- metadata_ordered %>% dplyr::filter(class %in% c("first_line", "second_line")) %>% pull(description)
#poi <- stats_res_tib_poi %>% dplyr::filter(feature %in% p_all) %>% pull(program) %>% unique()
