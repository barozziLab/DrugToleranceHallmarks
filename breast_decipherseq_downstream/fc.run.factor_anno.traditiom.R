
#~~~~~~~~~~~~~~~~~~~~~~
#Prepare TRADITIOM data
#######################

trad_home <- "TRADITIOM/"

#bulk-RNA-seq: VST normalised counts
inF <- paste0(trad_home, "bulk_RNAseq/vsd.counts.txt")
bulk_vst_counts <- read.delim(inF, sep = "\t", header = TRUE)
bulk_vst_counts_tib <- bulk_vst_counts %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble()

#scRNA-seq
inF <- paste0(trad_home, "/singlets.combinedAI.MNN0.4.rds")
so <- readRDS(inF)
#orig.ident as ordered factor
so@meta.data$orig.ident <- factor(so@meta.data$orig.ident, 
                                  levels = c("T0_1", "T0_2", 
                                             "AI1m1", "AI2m1", "AI1m2", "AI2m2", 
                                             "AIaw1", "AIaw2", "AIaw3", "AIaw4"))

#note: same top number of genes used for the enrichment analyses
#i.e. top_genes_n


#~~~~~~~~~~~~~~~~~~~~~
#ssGSEA (bulk-RNA-seq)
######################

#https://www.rdocumentation.org/packages/corto/versions/1.1.0/topics/ssgsea

#merge programs markers with PA signatures
sets_oi <- mrk_list
names(sets_oi) <- paste0("Program_", names(sets_oi))
pa_gsets <- m_t2g$PA$gs_name %>% unique
for (pa_gset in pa_gsets) {
  sets_oi[[pa_gset]] <- m_t2g$PA %>% 
    dplyr::filter(gs_name == pa_gset) %>% 
    pull(gene_symbol)
}

#run ssGSEA
nesmat <- ssgsea(bulk_vst_counts, sets_oi)
nesdf <- nesmat %>% as.data.frame()
colnames(nesdf) <- colnames(bulk_vst_counts)

#heat map
col_anno <- c(rep("POT", 3), rep("UT", 7), 
              rep("Dormancy_TAM", 4), rep("Dormancy_AI", 4), rep("Dormancy_TAM", 4), rep("Dormancy_AI", 6),
              rep("Awakening_TAM", 8), rep("Awakening_AI", 5),
              rep("TEP_TAM", 6), rep("TEP_AI", 4))
col_anno <- factor(col_anno, levels = unique(col_anno))
col_fun_scaled_nes <- colorRamp2(c(-10, 0, 10), colors = c("navy", "white", "red"))
set.seed(123)
hm <- Heatmap(nesdf, 
              col = col_fun_scaled_nes, 
              cluster_rows = TRUE, 
              clustering_distance_rows = "spearman", 
              cluster_columns = FALSE, 
              column_split = col_anno,
              name = "NES", 
              heatmap_legend_param = list(title = "NES"), 
              row_names_max_width = unit(10, "cm"), 
              border = TRUE)
out_path <- paste0(out_folder_traditiom, "/bulkRNAseq.ssGSEA.heatmap.pdf")
pdf(out_path, width = 14, height = 9)
draw(hm)
dev.off()

#save in long format
nestib_long <- nesdf %>% 
  rownames_to_column(var = "gs_name") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c("gs_name"), names_to = "sample", values_to = "NES")
#adding annotation
anno_tib <- tibble(sample = colnames(nesdf), annotation = col_anno)
nestib_long <- nestib_long %>% 
  left_join(anno_tib, by = "sample") %>% 
  select(gs_name, sample, annotation, NES)
out_path <- paste0(out_folder_traditiom, "/bulkRNAseq.ssGSEA.txt") 
write_tsv(x = nestib_long, file = out_path)

#box plots
bplot <- ggplot(nestib_long, aes(x=annotation, y=NES)) + 
  geom_boxplot() +
  coord_flip() +
  geom_hline(yintercept = 0) +
  facet_wrap(~gs_name, scales = "free") + 
  theme_bw()
out_path <- paste0(out_folder_traditiom, "/bulkRNAseq.ssGSEA.boxplots.pdf")
pdf(out_path, width = 15, height = 12)
plot(bplot)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~
#moduleScores (scRNA-seq)
#########################

gsets_oi <- names(sets_oi)

for (gset_id in gsets_oi) {
  foi <- list(c(sets_oi[[gset_id]]))
  so <- AddModuleScore(object = so, features = foi, ctrl = 100, name = gset_id)
  w <- colnames(so@meta.data) == paste0(gset_id, 1)
  colnames(so@meta.data)[w] <- paste0("score.", gset_id)
}

#violin plots
plt_vln <- list()
for (gset_id in gsets_oi) {
  gset_name <- paste0("score.", gset_id)
  plt_vln[[gset_id]] <- VlnPlot(so, features = gset_name, group.by = "orig.ident", pt.size = 0) + 
    ggtitle(gset_id) +
    stat_summary(fun.y = median, geom='point', size = 10, colour = "black", shape = 95)
}
out_path <- paste0(out_folder_traditiom, "/scRNAseq.moduleScore.vioplots.pdf")
pdf(out_path, width = 7, height = 5)
for (gset_id in gsets_oi) {
  plot(plt_vln[[gset_id]])
}
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#scRNA-seq: reconstruct lineage information
###########################################

#by lineage (barcode)

bc_table <- so@assays$Custom@data %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'barcode') %>% 
  as_tibble()

bc_table_long <- bc_table %>% 
  pivot_longer(cols = -c("barcode")) %>% 
  dplyr::filter(value > 0)

#exclude barcode that were counted only once in a cell, then compute stats per cell

bc_table_long_gr <- bc_table_long %>% 
  dplyr::filter(value > 1) %>%
  group_by(name) %>% 
  summarise(n = n(), sum = sum(value), max = max(value)) %>%
  mutate(max_frac = max / sum)

#retain only cells with max_frac >= 0.6 & sum >= 5
#some cells with multiple barcodes are retained, but they have to be quite clean (60% dominated)

bc_table_filter <- bc_table_long_gr %>% 
  filter(max_frac >= 0.6 & sum >= 5)

#consolidate cell <-> best barcode pairs
bc_table_max <- bc_table_long %>% 
  group_by(name) %>%
  filter(value == max(value)) %>%
  ungroup()
w <- bc_table_max$name %in% bc_table_filter$name
bc_table_max <- bc_table_max[w,]

#new so object (subset based on valid cells)
so_subset <- subset(x = so, sampleCell %in% bc_table_max$name)

#add info to the metadata, about the barcode, and the class (winning or other)
bc_winners <- c("bc-5122096", "bc-7704123", "bc-2127641", "bc-6864523")

metadata <- so_subset@meta.data %>% 
  as_tibble() %>% 
  left_join(bc_table_max, by=c("sampleCell" = "name")) %>%
  dplyr::rename(barcode_readcount = value) %>%
  mutate(barcode_winner = barcode %in% bc_winners)

so_subset@meta.data <- metadata %>% column_to_rownames("sampleCell")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#moduleScores (scRNA-seq; lineage)
##################################

#prepare data

w <- grepl("score", colnames(so_subset@meta.data))
cnames <- c("orig.ident", "barcode", "barcode_winner", colnames(so_subset@meta.data)[w])

bc_data_sig <- so_subset@meta.data[,cnames] %>% 
  rownames_to_column(var = "cell") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c("cell", "orig.ident", "barcode", "barcode_winner")) %>%
  mutate(barcode_class = ifelse(barcode_winner, barcode, "other"))

bc_data_sig$name <- gsub("score\\.", "", bc_data_sig$name)

#summarise
bc_data_sig_stats <- bc_data_sig %>% 
  group_by(barcode_class, name, orig.ident) %>% 
  summarise(value_mean = mean(value)) %>% 
  ungroup() %>%
  mutate(value_mean_sign = sign(value_mean)) 

#for the awakening, only include the cells with the winning barcode, specific for each awakening
#bc-2127641 -> AIaw3
bc_data_sig_stats <- bc_data_sig_stats %>% 
  dplyr::filter(!(barcode_class == "bc-2127641" & orig.ident %in% c("AIaw1", "AIaw2", "AIaw4")))
#bc-5122096 -> AIaw4
bc_data_sig_stats <- bc_data_sig_stats %>% 
  dplyr::filter(!(barcode_class == "bc-5122096" & orig.ident %in% c("AIaw1", "AIaw2", "AIaw3")))
#bc-6864523 -> AIaw2
bc_data_sig_stats <- bc_data_sig_stats %>% 
  dplyr::filter(!(barcode_class == "bc-6864523" & orig.ident %in% c("AIaw1", "AIaw3", "AIaw4")))
#bc-7704123 -> AIaw1
bc_data_sig_stats <- bc_data_sig_stats %>% 
  dplyr::filter(!(barcode_class == "bc-7704123" & orig.ident %in% c("AIaw2", "AIaw3", "AIaw4")))

#plots

programs <- bc_data_sig_stats$name %>% unique()
plts <- list()
for (program in programs) {
  plts[[program]] <- ggplot(bc_data_sig_stats %>% dplyr::filter(name == program), aes(x=orig.ident, y=value_mean, group=barcode_class)) +
    geom_line(aes(linetype=barcode_class)) +
    geom_point(aes(shape=barcode_class)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylim(-0.1, 1.4) +
    xlab("Identity") + 
    ylab("Mean Program Usage") +
    ggtitle(program)
}

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.lineplots.pdf")
pdf(out_path, width = 7, height = 5)
for (program in programs) {
  plot(plts[[program]])
}
dev.off()

plts_bar <- list()
for (program in programs) {
  plts_bar[[program]] <- ggplot(data=bc_data_sig_stats %>% dplyr::filter(name == program), aes(x=orig.ident, y=value_mean, fill=barcode_class)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "grey")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylim(-0.1, 1.4) +
    xlab("Identity") + 
    ylab("Mean Program Usage") +
    ggtitle(program)
}

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.barplots.pdf")
pdf(out_path, width = 7, height = 5)
for (program in programs) {
  plot(plts_bar[[program]])
}
dev.off()

#corresponding box plots

#for the awakening, only include the cells with the winning barcode, specific for each awakening
#bc-2127641 -> AIaw3
bc_data_sig_filt <- bc_data_sig %>% 
  dplyr::filter(!(barcode_class == "bc-2127641" & orig.ident %in% c("AIaw1", "AIaw2", "AIaw4")))
#bc-5122096 -> AIaw4
bc_data_sig_filt <- bc_data_sig_filt %>% 
  dplyr::filter(!(barcode_class == "bc-5122096" & orig.ident %in% c("AIaw1", "AIaw2", "AIaw3")))
#bc-6864523 -> AIaw2
bc_data_sig_filt <- bc_data_sig_filt %>% 
  dplyr::filter(!(barcode_class == "bc-6864523" & orig.ident %in% c("AIaw1", "AIaw3", "AIaw4")))
#bc-7704123 -> AIaw1
bc_data_sig_filt <- bc_data_sig_filt %>% 
  dplyr::filter(!(barcode_class == "bc-7704123" & orig.ident %in% c("AIaw2", "AIaw3", "AIaw4")))

plts_box <- list()
for (program in programs) {
  plts_box[[program]] <- bc_data_sig_filt %>% 
    dplyr::filter(name == program) %>% 
    ggplot(aes(x = orig.ident, y = value, fill = barcode_class)) + 
    geom_boxplot() +
    scale_fill_manual(values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "grey")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylim(-0.1, 1.4) +
    xlab("Identity") + 
    ylab("Program Usage") +
    ggtitle(program)
}

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.boxplots.pdf")
pdf(out_path, width = 7, height = 5)
for (program in programs) {
  plot(plts_box[[program]])
}
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~
#scRNA-seq; lineage -> ML
#########################

##Are there features correlated with successful awakenings?

ml_meta <- bc_data_sig %>% 
  dplyr::filter(grepl("AIaw", orig.ident)) %>% 
  select(cell, orig.ident, barcode_class, barcode_winner) %>% 
  unique()

ml_data <- bc_data_sig %>% 
  dplyr::filter(grepl("AIaw", orig.ident)) %>% 
  select(cell, name, value) %>% pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames('cell') %>%
  as.matrix()

#Train/Test on each AIaw separately, and save all results

#ident_ids <- c("AIaw1", "AIaw2", "AIaw3", "AIaw4")
#bc_wins <- c("bc-7704123", "bc-6864523", "bc-2127641", "bc-5122096")
#
#AIaw2 has too few negative examples, excluded
ident_ids <- c("AIaw1", "AIaw3", "AIaw4")
bc_wins <- c("bc-7704123", "bc-2127641", "bc-5122096")

ml_results <- list()
ml_results[["perf"]] <- list()
ml_results[["imp"]] <- list()

for (i in 1:length(ident_ids)) {
  
  ident_id <- ident_ids[i]
  bc_win <- bc_wins[i]
  
  #barcode_winner becomes y
  ml_meta_i <- ml_meta %>% 
    dplyr::filter(orig.ident == ident_id) %>%
    mutate(barcode_winner = ifelse(barcode_class == bc_win, TRUE, FALSE)) 
  
  ml_data_i <- ml_data %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell") %>% 
    as_tibble() %>% 
    left_join(ml_meta_i %>% select(cell, barcode_winner), by = "cell") %>%
    dplyr::rename(y = barcode_winner) %>%
    mutate(y = as.factor(y)) %>%
    select(-cell)
  
  #force manual subsampling
  set.seed(123)
  w_pos <- sample(which(ml_data_i$y == TRUE), 100)
  set.seed(123)
  w_neg <- sample(which(ml_data_i$y == FALSE), 100)
  ml_data_i_ss <- ml_data_i[c(w_pos, w_neg), ]
  
  set.seed(123)
  ml_split <- initial_split(ml_data_i_ss)
  ml_train <- training(ml_split)
  ml_test <- testing(ml_split)
  
  ml_recipe <- training(ml_split) %>%
    recipe(y ~.) %>%
    step_center(all_predictors(), -all_outcomes()) %>%
    step_scale(all_predictors(), -all_outcomes()) %>%
    prep()
  
  ## RF (simple, no tuning)
  #https://rviews.rstudio.com/2019/06/19/a-gentle-intro-to-tidymodels/
  
  ml_testing <- ml_recipe %>%
    bake(testing(ml_split)) 
  
  ml_training <- juice(ml_recipe)
  
  ml_ranger <- rand_forest(trees = 100, mode = "classification") %>%
    set_engine("ranger", importance = "permutation") %>%
    fit(y ~ ., data = ml_training)
  
  my_metrics <- metric_set(accuracy, sens, spec, kap)
  
  rf_confusion_mat <- ml_ranger %>%
    predict(ml_testing) %>%
    bind_cols(ml_testing) %>%
    my_metrics(truth = y, estimate = .pred_class)
  ml_results[["perf"]][[bc_win]] <- rf_confusion_mat
  
  #plt <- ml_ranger %>%
  #  vip(geom = "col", num_features = 40) + theme_bw()
  
  rf_feature_importance <- vip:::vi(ml_ranger)
  ml_results[["imp"]][[bc_win]] <- rf_feature_importance

}

#consolidate the results

res_perf <- c()
res_imp <- c()

for (i in 1:length(ident_ids)) {
  
  ident_id <- ident_ids[i]
  bc_win <- bc_wins[i]
  
  res_perf <- rbind(res_perf, ml_results[["perf"]][[bc_win]] %>% 
                      mutate(ident = ident_id, barcode = bc_win))
  res_imp <- rbind(res_imp, ml_results[["imp"]][[bc_win]] %>% 
                     mutate(ident = ident_id, barcode = bc_win) %>% 
                     arrange(-Importance) %>% 
                     mutate(Rank = 1:nrow(ml_results[["imp"]][[bc_win]])))

}

#sort by average importance
lvls <- res_imp %>% 
  group_by(Variable) %>% 
  summarise(avg = mean(Importance)) %>% 
  arrange(-avg) %>% 
  pull(Variable) %>%
  rev()
res_imp <- res_imp %>%
  mutate(Variable = factor(Variable, levels = lvls)) %>%
  mutate(Importance_sign = as.factor(sign(Importance)))

plt_imp <- ggplot(res_imp, aes(x = Variable, y = ident)) + 
  geom_point(aes(size = Importance, fill = Rank), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "red", high = "white") + 
  theme(legend.key.size = unit(0.4, 'cm')) +
  coord_flip()

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML.importance.pdf")
pdf(out_path, width = 3.5, height = 6.5)
plot(plt_imp)
dev.off()

plt_perf <- ggplot(res_perf, aes(x = .metric, y = .estimate)) + 
  geom_boxplot() +
  geom_jitter(width = 0.05) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0.4, 1) +
  xlab("Metric") +
  ylab("Estimate") +
  ggtitle("Performances")

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML.performance.pdf")
pdf(out_path, width = 2, height = 3)
plot(plt_perf)
dev.off()

#load the results of functional enrichment analyses
#extract and plot the top 3 programs

gsea_res <- read_tsv("fc_results/W.enrichr.results.txt")

top_programs <- res_imp$Variable %>% head(3)
top_programs <- gsub("Program_", "", top_programs)

top_terms <- 30
ylim_manual <- c(0, 12)

plt_gsea <- list()

for (i in 1:length(top_programs)) {

  plt_res <- gsea_res %>% 
    dplyr::filter(Gene_set == top_programs[i]) %>%
    arrange(p.adjust) %>% 
    head(30)

  #sort the results (the IDs) by p.adjust
  lvls <- plt_res %>% arrange(p.adjust) %>% pull(ID) %>% rev()
  plt_res <- plt_res %>% mutate(ID = factor(ID, levels = lvls))
  
  plt_gsea[[top_programs[i]]] <- ggplot(data = plt_res, aes(x=ID, y=-log10(p.adjust), fill=Count)) + 
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05)) +
    coord_flip() +
    theme_bw() +
    xlab("") +
    ggtitle(top_programs[i]) +
    ylim(ylim_manual)

}
out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML.importance.GSEA.top_30.pdf")
pdf(out_path, width = 6, height = 6)
for (i in 1:length(top_programs)) {
  plot(plt_gsea[[top_programs[i]]])
}
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#scRNA-seq; lineage -> ML all --> winner barcodes only
######################################################

##Are there features correlated with successful barcodes, even before awakenings?
##ML on T0, AI, Aw, separately

ml_meta <- bc_data_sig %>% 
  select(cell, orig.ident, barcode_class, barcode_winner) %>% 
  unique()

ml_data <- bc_data_sig %>% 
  select(cell, name, value) %>% pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames('cell') %>%
  as.matrix()

#Train/Test on each sample separately, and save all results
i_ids <- c("T0_1", "T0_2", "AI1m1", "AI2m1", "AI1m2", "AI2m2", "AIaw1", "AIaw2", "AIaw3", "AIaw4")
ident_ids <- c()
b_ids <- c("bc-7704123", "bc-6864523", "bc-2127641", "bc-5122096")
bc_wins <- c()
for (i_id in i_ids) {
  for (b_id in b_ids) {
    ident_ids <- c(ident_ids, i_id)
    bc_wins <- c(bc_wins, b_id)
  }
}

ml_results <- list()
ml_results[["perf"]] <- list()
ml_results[["imp"]] <- list()
ml_results[["pos_class"]] <- list() # contains the ids of the cells in the test sets that are classified as positive (can be used to fish out outliers from T0 and AI)

Ks <- 1:10

for (i in 1:length(ident_ids)) {
  
  ident_id <- ident_ids[i]
  bc_win <- bc_wins[i]
    
  #barcode_winner becomes y
  ml_meta_i <- ml_meta %>% 
    dplyr::filter(orig.ident == ident_id) %>%
    mutate(barcode_winner = ifelse(barcode_class == bc_win, TRUE, FALSE))
    
  ml_data_i <- ml_data %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell") %>% 
    as_tibble() %>% 
    inner_join(ml_meta_i %>% select(cell, barcode_winner), by = "cell") %>%
    dplyr::rename(y = barcode_winner) %>%
    mutate(y = as.factor(y))
  
  pos_examples <- sum(ml_data_i$y == TRUE)
  neg_examples <- sum(ml_data_i$y == FALSE)
  
  for (k in Ks) {
    
    combined_id <- paste0(ident_id, "__", bc_win, "__", k)
    
    if (pos_examples >= 100 & neg_examples >= 100) {
      
      #force manual subsampling
      set.seed(k)
      w_pos <- sample(which(ml_data_i$y == TRUE), 100)
      set.seed(k)
      w_neg <- sample(which(ml_data_i$y == FALSE), 100)
      ml_data_i_ss <- ml_data_i[c(w_pos, w_neg), ]
      ml_data_i_ss <- ml_data_i_ss %>% column_to_rownames("cell")
      
      set.seed(123)
      ml_split <- initial_split(ml_data_i_ss)
      ml_train <- training(ml_split)
      ml_test <- testing(ml_split)
      
      ml_recipe <- training(ml_split) %>%
        recipe(y ~.) %>%
        step_center(all_predictors(), -all_outcomes()) %>%
        step_scale(all_predictors(), -all_outcomes()) %>%
        prep()
      
      ## RF (simple, no tuning)
      #https://rviews.rstudio.com/2019/06/19/a-gentle-intro-to-tidymodels/
      
      ml_testing <- ml_recipe %>%
        bake(testing(ml_split))
      
      ml_training <- juice(ml_recipe)
      
      ml_ranger <- rand_forest(trees = 100, mode = "classification") %>%
        set_engine("ranger", importance = "permutation") %>%
        fit(y ~ ., data = ml_training)
      
      my_metrics <- metric_set(accuracy, sens, spec, kap)
      
      rf_confusion_mat <- ml_ranger %>%
        predict(ml_testing) %>%
        bind_cols(ml_testing) %>%
        my_metrics(truth = y, estimate = .pred_class)
      ml_results[["perf"]][[combined_id]] <- rf_confusion_mat
      
      #plt <- ml_ranger %>%
      #  vip(geom = "col", num_features = 40) + theme_bw()
      
      rf_feature_importance <- vip:::vi(ml_ranger)
      ml_results[["imp"]][[combined_id]] <- rf_feature_importance
      
      #save ids of cells in the test sets with a score >= 0.9
      
      score_t <- 0.9
      
      ml_testing_prob <- ml_ranger %>% 
        predict(ml_testing, type="prob") %>%
        mutate(consider = .pred_TRUE >= score_t)
      consider_i <- which(ml_testing_prob$consider)
      pos_class_ids <- rownames(ml_test[consider_i,])
      
      ml_results[["pos_class"]][[combined_id]] <- c(ml_results[["pos_class"]][[combined_id]], pos_class_ids)
    
    }
    
  }
  
}

#consolidate the results

res_perf <- c()
res_imp <- c()
res_posCls <- c()

for (combined_id in names(ml_results[["perf"]])) {
  
  ident_id <- strsplit(combined_id, "__")[[1]][1]
  bc_id <- strsplit(combined_id, "__")[[1]][2]
  sample_k <- strsplit(combined_id, "__")[[1]][3]
  
  res_perf <- rbind(res_perf, ml_results[["perf"]][[combined_id]] %>% 
                      mutate(ident = ident_id, barcode = bc_id, k = sample_k))
  
  res_imp <- rbind(res_imp, ml_results[["imp"]][[combined_id]] %>% 
                     mutate(ident = ident_id, barcode = bc_id, k = sample_k, id_combined = combined_id) %>% 
                     arrange(-Importance) %>% 
                     mutate(Rank = 1:nrow(ml_results[["imp"]][[combined_id]])))
  
  if (length(ml_results[["pos_class"]][[combined_id]]) > 0) {
    
    res_posCls <- rbind(res_posCls,  
                        tibble(ident = ident_id, barcode = bc_id, k = sample_k, id_combined = combined_id, 
                               cell_id = ml_results[["pos_class"]][[combined_id]]))
    
  }
  
}

res_perf <- res_perf %>% mutate(ident = factor(ident, levels = i_ids))
res_perf <- res_perf %>% mutate(barcode = factor(barcode, levels = b_ids))
res_perf <- res_perf %>% mutate(k = as.integer(k))

res_imp <- res_imp %>% mutate(ident = factor(ident, levels = i_ids))
res_imp <- res_imp %>% mutate(barcode = factor(barcode, levels = b_ids))
res_imp <- res_imp %>% mutate(k = as.integer(k))

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.importance.full.txt")
write_tsv(x = res_imp, file = out_path)

#classified as positive

#save
out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.testing_pos_class.txt")
write_tsv(x = res_posCls, file = out_path)

#statistics
res_posCls_stats <- res_posCls %>% 
  group_by(ident, barcode, k) %>% 
  summarise(n = n()) %>% 
  mutate(fraction_approx = n / 25) %>% 
  ungroup() %>% 
  mutate(ident = factor(ident, levels = i_ids))

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.testing_pos_class.stats.txt")
write_tsv(x = res_posCls_stats, file = out_path)

plt_posCls_stats <- res_posCls_stats %>%
  ggplot(aes(x = ident, y = fraction_approx)) + 
  geom_boxplot() +
  theme_bw() +
  ylab("Fraction (approx.)") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.testing_pos_class.pdf")
pdf(out_path, width = 2.5, height = 2.5)
plot(plt_posCls_stats)
dev.off()

## Performances

plt_perf <- ggplot(res_perf %>% dplyr::filter(.metric != "kap"), aes(x = ident, y = .estimate)) + 
  geom_boxplot() +
  facet_grid(barcode ~ .metric) +
  #geom_jitter(width = 0.05) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0.2, 1) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 0.6, col = "red") +
  geom_hline(yintercept = 0.7, col = "orange") +
  geom_hline(yintercept = 0.8, col = "darkgreen") +
  xlab("Sample") +
  ylab("Estimate") +
  ggtitle("Performances")

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.performance.pdf")
pdf(out_path, width = 5, height = 6)
plot(plt_perf)
dev.off()

## Heat map of importance

res_imp_mat <- res_imp %>% 
  select(Variable, id_combined, Importance) %>% 
  pivot_wider(names_from = Variable, values_from = Importance) %>% 
  column_to_rownames("id_combined")

#filter on features
importance_t <- 0.02
col_w <- apply(res_imp_mat, 2, max) >= importance_t
res_imp_mat <- res_imp_mat[,col_w]

row_anno_condition <- sapply(rownames(res_imp_mat), function(x){gsub("[1-4]$|_[1-4]$", "", strsplit(x, split = "__")[[1]][1])}) %>% as.character()
row_anno_sample <- sapply(rownames(res_imp_mat), function(x){gsub("T0_|AIm|AI1m|AI2m|AIaw", "", strsplit(x, split = "__")[[1]][1])}) %>% as.character()

row_anno_condition_lvls <- c("T0", "AI1m", "AI2m", "AIaw")
row_anno_condition <- factor(row_anno_condition, levels = row_anno_condition_lvls)
row_anno_condition_cols <- RColorBrewer::brewer.pal(name = "Accent", n = 4)
names(row_anno_condition_cols) <- row_anno_condition_lvls
row_anno_sample_lvls <- c("1", "2", "3", "4")
row_anno_sample <- factor(row_anno_sample, levels = row_anno_sample_lvls)
row_anno_sample_cols <- RColorBrewer::brewer.pal(name = "Blues", n = 4)
names(row_anno_sample_cols) <- row_anno_sample_lvls

ha <- rowAnnotation(condition = row_anno_condition, 
                    sample_n = row_anno_sample, 
                    col = list(condition = row_anno_condition_cols, 
                               sample_n = row_anno_sample_cols))

set.seed(123)
hm <- Heatmap(as.matrix(res_imp_mat), 
              cluster_rows = TRUE, 
              cluster_columns = TRUE, 
              right_annotation = ha,
              name = "Importance", 
              heatmap_legend_param = list(title = "Importance"), 
              row_names_max_width = unit(10, "cm"), 
              row_names_gp = gpar(fontsize = 2), 
              border = TRUE)

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.importance.heatmap.pdf")
pdf(out_path, width = 8, height = 9)
draw(hm)
dev.off()

## Bubble plot of mean importance per feature per condition

res_imp_summary <- res_imp %>% 
  group_by(Variable, ident) %>% 
  summarise(Importance = mean(Importance), Rank = mean(Rank)) %>% 
  ungroup()

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.importance.summary.txt")
write_tsv(x = res_imp_summary, file = out_path)

#remove programs with no importance(s) over a certain threshold (not that low)
imp_t <- 0.01
res_features <- res_imp_summary %>% 
  group_by(Variable) %>% 
  summarise(Importance_max = max(Importance)) %>% 
  dplyr::filter(Importance_max >= imp_t) %>%
  pull(Variable)
res_imp_summary_plt <- res_imp_summary %>% 
  dplyr::filter(Variable %in% res_features)

#set all importance(s) below a certain (low) threshold to NAs (to avoid displaying them)
imp_t <- 0.002
res_imp_summary_plt <- res_imp_summary_plt %>% 
  mutate(Importance = ifelse(Importance <= imp_t, NA, Importance))

plt_imp <- ggplot(res_imp_summary_plt, aes(x = Variable, y = ident)) + 
  geom_point(aes(size = Importance, fill = Rank), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "red", high = "white") + 
  theme(legend.key.size = unit(0.4, 'cm')) +
  coord_flip()

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.importance.average.pdf")
pdf(out_path, width = 4.5, height = 3.5)
plot(plt_imp)
dev.off()

#load the results of functional enrichment analyses
#extract and plot selected programs

gsea_res <- read_tsv("fc_results/W.enrichr.results.txt")

selected_programs <- c("Program_25", "Program_18", "Program_14", "Program_12", "Program_10")
selected_programs <- gsub("Program_", "", selected_programs)

top_terms <- 30
ylim_manual <- c(0, 12)

plt_gsea <- list()

for (i in 1:length(selected_programs)) {
  
  plt_res <- gsea_res %>% 
    dplyr::filter(Gene_set == selected_programs[i]) %>%
    arrange(p.adjust) %>% 
    head(30)
  
  #sort the results (the IDs) by p.adjust
  lvls <- plt_res %>% arrange(p.adjust) %>% pull(ID) %>% rev()
  plt_res <- plt_res %>% mutate(ID = factor(ID, levels = lvls))
  
  #saturate p.adjust
  plt_res_sat <- plt_res %>% mutate(p.adjust = ifelse(p.adjust <= 1e-12, 1e-12, p.adjust))
  
  plt_gsea[[selected_programs[i]]] <- ggplot(data = plt_res_sat, aes(x=ID, y=-log10(p.adjust), fill=Count)) + 
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05)) +
    coord_flip() +
    theme_bw() +
    xlab("") +
    ggtitle(selected_programs[i]) +
    ylim(ylim_manual)
  
}
out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.importance.average.GSEA.top_30.individual.pdf")
pdf(out_path, width = 6, height = 6)
for (i in 1:length(selected_programs)) {
  plot(plt_gsea[[selected_programs[i]]])
}
dev.off()

#bubble plot of the same programs, only top 10 terms, excluding the cell-type signatures

ids_merge <- c()

for (i in 1:length(selected_programs)) {
  
  ids <- gsea_res %>% 
    dplyr::filter(Gene_set == selected_programs[i]) %>%
    dplyr::filter(!grepl("cell-state_", ID)) %>%
    arrange(p.adjust) %>% 
    pull(ID) %>%
    head(10)
  
  #save the top ids
  ids_merge <- c(ids_merge, ids)
  
}

selected_programs_sorted <- c("14", "12", "10", "25", "18")

gsea_res_input <- gsea_res %>% 
  dplyr::filter(Gene_set %in% selected_programs & ID %in% ids_merge) %>% 
  mutate(Gene_set = factor(Gene_set, levels = selected_programs_sorted)) %>%
  mutate(p.adjust = ifelse(p.adjust <= 1e-10, 1e-10, p.adjust))

#matrix for hclust (for re-ordering)
gsea_res_mat <- gsea_res_input %>% 
  select(Gene_set, ID, p.adjust) %>% 
  mutate(p.adjust = -log10(p.adjust)) %>% 
  pivot_wider(names_from = Gene_set, values_from = p.adjust) %>%
  column_to_rownames("ID")
gsea_res_mat[is.na(gsea_res_mat)] <- 0

hc_programs <- hclust(dist(gsea_res_mat))

lvls_sorted <- rownames(gsea_res_mat)[hc_programs$order]

gsea_res_input <- gsea_res_input %>%
  mutate(ID = factor(ID, levels = rev(lvls_sorted)))

plt_imp <- gsea_res_input %>%
  ggplot(aes(x = ID, y = Gene_set)) + 
  geom_point(aes(size = Count, fill = -log10(p.adjust)), alpha = 0.75, shape = 21) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_gradient(low = "white", high = "darkblue") + 
  theme(legend.key.size = unit(0.4, 'cm')) +
  coord_flip()

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_YES.importance.average.GSEA.top_10.together.pdf")
pdf(out_path, width = 6.2, height = 6)
plot(plt_imp)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#scRNA-seq; lineage -> ML all --> no winners
############################################

## Expected distribution from all the other barcodes: try all barcodes with enough cells, excluding the winner

exp_ml_meta <- bc_data_sig %>% 
  select(cell, orig.ident, barcode) %>% 
  unique()

exp_ml_data <- bc_data_sig %>% 
  select(cell, name, value) %>% pivot_wider(names_from = "name", values_from = "value") %>%
  column_to_rownames('cell') %>%
  as.matrix()

#Train/Test on each sample separately, and save all results
i_ids <- c("T0_1", "T0_2", "AI1m1", "AI2m1", "AI1m2", "AI2m2")
ident_ids <- c()
b_ids <- exp_ml_meta$barcode %>% unique()
#exclude winnders
b_ids <- b_ids[!b_ids %in% c("bc-7704123", "bc-6864523", "bc-2127641", "bc-5122096")]
bc_wins <- c()
for (i_id in i_ids) {
  for (b_id in b_ids) {
    ident_ids <- c(ident_ids, i_id)
    bc_wins <- c(bc_wins, b_id)
  }
}

exp_ml_results <- list()
exp_ml_results[["perf"]] <- list()
exp_ml_results[["imp"]] <- list()
exp_ml_results[["pos_class"]] <- list() # contains the ids of the cells in the test sets that are classified as positive (can be used to fish out outliers from T0 and AI)

Ks <- 1:10

for (i in 1:length(ident_ids)) {
  
  ident_id <- ident_ids[i]
  bc_win <- bc_wins[i]
  
  #barcode_winner becomes y
  ml_meta_i <- exp_ml_meta %>% 
    dplyr::filter(orig.ident == ident_id) %>%
    mutate(barcode_winner = ifelse(barcode == bc_win, TRUE, FALSE))
  
  ml_data_i <- exp_ml_data %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell") %>% 
    as_tibble() %>% 
    inner_join(ml_meta_i %>% select(cell, barcode_winner), by = "cell") %>%
    dplyr::rename(y = barcode_winner) %>%
    mutate(y = as.factor(y))
  
  pos_examples <- sum(ml_data_i$y == TRUE)
  neg_examples <- sum(ml_data_i$y == FALSE)
  
  for (k in Ks) {
    
    combined_id <- paste0(ident_id, "__", bc_win, "__", k)
    
    if (pos_examples >= 100 & neg_examples >= 100) {
      
      #force manual subsampling
      set.seed(k)
      w_pos <- sample(which(ml_data_i$y == TRUE), 100)
      set.seed(k)
      w_neg <- sample(which(ml_data_i$y == FALSE), 100)
      ml_data_i_ss <- ml_data_i[c(w_pos, w_neg), ]
      ml_data_i_ss <- ml_data_i_ss %>% column_to_rownames("cell")
      
      set.seed(123)
      ml_split <- initial_split(ml_data_i_ss)
      ml_train <- training(ml_split)
      ml_test <- testing(ml_split)
      
      ml_recipe <- training(ml_split) %>%
        recipe(y ~.) %>%
        step_center(all_predictors(), -all_outcomes()) %>%
        step_scale(all_predictors(), -all_outcomes()) %>%
        prep()
      
      ## RF (simple, no tuning)
      #https://rviews.rstudio.com/2019/06/19/a-gentle-intro-to-tidymodels/
      
      ml_testing <- ml_recipe %>%
        bake(testing(ml_split))
      
      ml_training <- juice(ml_recipe)
      
      ml_ranger <- rand_forest(trees = 100, mode = "classification") %>%
        set_engine("ranger", importance = "permutation") %>%
        fit(y ~ ., data = ml_training)
      
      my_metrics <- metric_set(accuracy, sens, spec, kap)
      
      rf_confusion_mat <- ml_ranger %>%
        predict(ml_testing) %>%
        bind_cols(ml_testing) %>%
        my_metrics(truth = y, estimate = .pred_class)
      exp_ml_results[["perf"]][[combined_id]] <- rf_confusion_mat
      
      #plt <- ml_ranger %>%
      #  vip(geom = "col", num_features = 40) + theme_bw()
      
      rf_feature_importance <- vip:::vi(ml_ranger)
      exp_ml_results[["imp"]][[combined_id]] <- rf_feature_importance
      
    }
    
  }
  
}

#consolidate the results

exp_res_perf <- c()
exp_res_imp <- c()

for (combined_id in names(exp_ml_results[["perf"]])) {
  
  ident_id <- strsplit(combined_id, "__")[[1]][1]
  bc_id <- strsplit(combined_id, "__")[[1]][2]
  sample_k <- strsplit(combined_id, "__")[[1]][3]
  
  exp_res_perf <- rbind(exp_res_perf, exp_ml_results[["perf"]][[combined_id]] %>% 
                        mutate(ident = ident_id, barcode = bc_id, k = sample_k))
  
  exp_res_imp <- rbind(exp_res_imp, exp_ml_results[["imp"]][[combined_id]] %>% 
                       mutate(ident = ident_id, barcode = bc_id, k = sample_k, id_combined = combined_id) %>% 
                       arrange(-Importance) %>% 
                       mutate(Rank = 1:nrow(exp_ml_results[["imp"]][[combined_id]])))
  
}

exp_res_perf <- exp_res_perf %>% mutate(ident = factor(ident, levels = i_ids))
exp_res_perf <- exp_res_perf %>% mutate(barcode = factor(barcode, levels = b_ids))
exp_res_perf <- exp_res_perf %>% mutate(k = as.integer(k))

exp_res_imp <- exp_res_imp %>% mutate(ident = factor(ident, levels = i_ids))
exp_res_imp <- exp_res_imp %>% mutate(barcode = factor(barcode, levels = b_ids))
exp_res_imp <- exp_res_imp %>% mutate(k = as.integer(k))

## Performances (non-winners)

plt_perf <- ggplot(exp_res_perf %>% dplyr::filter(.metric != "kap"), aes(x = ident, y = .estimate)) + 
  geom_boxplot() +
  facet_wrap(~.metric) +
  #geom_jitter(width = 0.05) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0.2, 1) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 0.6, col = "red") +
  geom_hline(yintercept = 0.7, col = "orange") +
  geom_hline(yintercept = 0.8, col = "darkgreen") +
  xlab("Sample") +
  ylab("Estimate") +
  ggtitle("Performances")

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.winners_NO.performance.pdf")
pdf(out_path, width = 5, height = 3)
plot(plt_perf)
dev.off()

##merge winners and non-winners

all_res_perf <- rbind(exp_res_perf %>% mutate(winner = "no"), res_perf %>% mutate(winner = "yes"))

plt_perf <- ggplot(all_res_perf %>% dplyr::filter(.metric != "kap"), aes(x = ident, y = .estimate, fill = winner)) + 
  geom_boxplot() +
  facet_wrap(~.metric) +
  #geom_jitter(width = 0.05) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylim(0.2, 1) +
  geom_hline(yintercept = 0.5) +
  geom_hline(yintercept = 0.6, col = "red") +
  geom_hline(yintercept = 0.7, col = "orange") +
  geom_hline(yintercept = 0.8, col = "darkgreen") +
  scale_fill_manual(values = c("grey", "lightblue")) +
  xlab("Sample") +
  ylab("Estimate") +
  ggtitle("Performances")

out_path <- paste0(out_folder_traditiom, "/scRNAseq.by_bc.moduleScore.ML_ALL.both.performance.pdf")
pdf(out_path, width = 7, height = 3)
plot(plt_perf)
dev.off()

