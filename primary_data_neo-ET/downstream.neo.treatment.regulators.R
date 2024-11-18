
source("downstream.libs.R")

outFolder <- "downstream.neo.treatment.results/"

ha_cols <- list(treatment = c("pre" = "grey", "post" = "purple"),
                response = c("Responder" = "darkgreen", "Non_responder" = "orange", "Stable" = "lightblue", "NA" = "lightgrey"))


#~~~~~~~~~~~~~~~~
#Data Preparation 
#################

dset_ens <- readRDS("sleuth.results/sleuth.results.expr_matrix.log2.rds")

#collapse to symbols

gene_map <- read_tsv("TranscriptToGene.txt")
colnames(gene_map) <- c("ensembl_gene", "ensembl_tx", "symbol", "ensembl_tx_v")

dset_ens <- dset_ens %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  dplyr::rename(ensembl_tx_v = rowname)

dset_ens_long <- dset_ens %>% 
  pivot_longer(!ensembl_tx_v, names_to = "sample", values_to = "expr") %>%
  right_join(gene_map, by = "ensembl_tx_v") %>%
  dplyr::filter(!is.na(sample))

dset_ens_long_avg <- dset_ens_long %>% 
  group_by(sample, symbol) %>% 
  summarise(expr_avg = mean(expr))

dset_tib <- dset_ens_long_avg %>% 
  pivot_wider(names_from = sample, values_from = expr_avg)

dset <- dset_tib %>% 
  select(-c("symbol")) %>% 
  as.data.frame()
rownames(dset) <- dset_tib %>% pull(symbol)

#add metadata for heat map annotation

meta_clinic <- read_tsv("metadata_clinical.txt")
meta_sample <- read_tsv("metadata_samples.txt")

meta_sample <- meta_sample %>% 
  rowwise %>% 
  mutate(id = paste(patient_id, "_", treatment, sep = "")) %>% 
  ungroup()

meta_full <- meta_sample %>% 
  left_join(meta_clinic, by = "patient_id") %>%
  dplyr::rename(response = `Responder/Non_responder/Stable`)

meta_full$response[is.na(meta_full$response)] <- "NA"

dset <- dset[, meta_full$id]

## also DEG call

degs <- read_tsv("sleuth.results/sleuth.results.txt")
degs_b_mat <- degs %>% 
  select(target_id, b) %>% 
  column_to_rownames("target_id")


#~~~~~~~~~~~~~
#DecoupleR TFs
##############

net_tfs <- get_collectri(organism='human', split_complexes=FALSE)

#sample activities
sample_acts_tfs <- run_ulm(mat=dset, 
                           net=net_tfs, 
                           .source='source', 
                           .target='target', 
                           .mor='mor', 
                           minsize = 10)

sample_acts_mat <- sample_acts_tfs %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

#activities directly from DEGs
#use b value from differential expression
contrast_acts <- run_ulm(mat=degs_b_mat, 
                         net=net_tfs, 
                         .source='source', 
                         .target='target', 
                         .mor='mor', minsize = 10)
#correct for multiple hypotheses testing
contrast_acts <- contrast_acts %>%
  mutate(p_value_BH = p.adjust(contrast_acts$p_value, method = "BH"))

#only significant ones, transpose and scale
w <- contrast_acts %>% 
  dplyr::filter(p_value_BH <= 0.05) %>%
  pull(source)

# Filter top TFs (only drawing from the significant ones) with either sign
n_tfs <- 25
f_contrast_acts <- contrast_acts %>%
  dplyr::filter(source %in% w) %>%
  mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)
f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)

#extract the matrix of this TFs for
sample_acts_mat_diff <- sample_acts_mat[,tfs]

#transpose and scale
sample_acts_mat_diff <- t(sample_acts_mat_diff)
cnames <- colnames(sample_acts_mat_diff)
sample_acts_mat_diff <- t(apply(sample_acts_mat_diff, 1, scale))
colnames(sample_acts_mat_diff) <- cnames

ha <- HeatmapAnnotation(treatment = meta_full[["treatment"]], 
                        response = meta_full[["response"]],
                        col = ha_cols)

HM <- Heatmap(sample_acts_mat_diff,
              heatmap_legend_param = list(title = "TF Activity (z-score)"),
              col = colorRamp2(c(-2,0,2), c("navy", "white", "red")),
              row_gap = unit(c(3), "mm"),
              row_names_gp = grid::gpar(fontsize = 8),
              column_names_gp = grid::gpar(fontsize = 8),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(3, "cm"),
              border = TRUE,
              top_annotation = ha,
              clustering_distance_rows = "pearson",
              clustering_distance_columns = "pearson")

outF <- paste(outFolder,"/neo.decoupleR.tfs.heatmap.pdf", sep = "")
pdf(outF, width=9, height=5)
draw(HM)
dev.off()

#same as box plots
sample_acts_filt <- sample_acts_mat[,tfs] %>% 
  as.data.frame %>% 
  rownames_to_column(var = "id") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c("id"), names_to = "pathway") %>% 
  rowwise() %>%
  mutate(patient = strsplit(id, split = "_")[[1]][1]) %>%
  mutate(treatment = strsplit(id, split = "_")[[1]][2])

bplot <- ggplot(data = sample_acts_filt, aes(x=treatment, y=value)) +
  facet_wrap(~pathway, ncol = 1) +
  geom_boxplot(outlier.size = 0) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.1) +
  coord_flip() +
  theme_bw()

outF <- paste(outFolder,"/neo.decoupleR.tfs.boxplots.pdf", sep = "")
pdf(outF, width=2, height=18)
plot(bplot)
dev.off()

#save a table of full results
outF <- paste(outFolder,"/neo.decoupleR.tfs.txt", sep = "")
write_tsv(x = contrast_acts, file = outF)


#~~~~~~~~~~~~~~~~~~~
#DecoupleR Signaling
####################

net_sig <- get_progeny(organism = 'human', top = 500)

#sample activities
sample_acts_sig <- run_ulm(mat=dset, 
                           net=net_sig, 
                           .source='source', 
                           .target='target', 
                           .mor='weight',
                           minsize = 5)

sample_acts_mat <- sample_acts_sig %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

#activities directly from DEGs
#use b value from differential expression
contrast_acts <- run_ulm(mat=degs_b_mat, 
                         net=net_sig, 
                         .source='source', 
                         .target='target', 
                         .mor='weight', 
                         minsize = 5)
#correct for multiple hypotheses testing
contrast_acts <- contrast_acts %>%
  mutate(p_value_BH = p.adjust(contrast_acts$p_value, method = "BH"))

#only significant ones, transpose and scale
w <- contrast_acts %>% 
  dplyr::filter(p_value_BH <= 0.05) %>%
  pull(source)
sample_acts_mat_scaled <- t(sample_acts_mat[,w])
sample_acts_mat_scaled <- t(apply(sample_acts_mat_scaled, 1, scale))
colnames(sample_acts_mat_scaled) <- rownames(sample_acts_mat)

ha <- HeatmapAnnotation(treatment = meta_full[["treatment"]], 
                        response = meta_full[["response"]],
                        col = ha_cols)

HM <- Heatmap(sample_acts_mat_scaled,
              heatmap_legend_param = list(title = "Sig. Pathway\n Activity (z-score)"),
              col = colorRamp2(c(-2,0,2), c("navy", "white", "red")),
              row_gap = unit(c(3), "mm"),
              row_names_gp = grid::gpar(fontsize = 8),
              column_names_gp = grid::gpar(fontsize = 8),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(3, "cm"),
              border = TRUE,
              top_annotation = ha)

outF <- paste(outFolder,"/neo.decoupleR.signaling.heatmap.pdf", sep = "")
pdf(outF, width=9, height=3.75)
draw(HM)
dev.off()

#same as box plots
sample_acts_filt <- sample_acts_mat[,w] %>% 
  as.data.frame %>% 
  rownames_to_column(var = "id") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c("id"), names_to = "pathway") %>% 
  rowwise() %>%
  mutate(patient = strsplit(id, split = "_")[[1]][1]) %>%
  mutate(treatment = strsplit(id, split = "_")[[1]][2])

bplot <- ggplot(data = sample_acts_filt, aes(x=treatment, y=value)) +
  facet_wrap(~pathway, ncol = 1) +
  geom_boxplot(outlier.size = 0) +
  geom_jitter(color="black", size=0.4, alpha=0.9, width = 0.1) +
  coord_flip() +
  ylim(-7.5, 7.5) +
  theme_bw()

outF <- paste(outFolder,"/neo.decoupleR.signaling.boxplots.pdf", sep = "")
pdf(outF, width=2, height=7)
plot(bplot)
dev.off()

#save a table of full results
outF <- paste(outFolder,"/neo.decoupleR.signaling.txt", sep = "")
write_tsv(x = contrast_acts, file = outF)


#~~~~~~~~~~~
#ML Modeling
############

#Signaling (no pre-filtering)
mat_sig <- sample_acts_sig %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>% t()

#TFs (pre-filtering on variance)
mat_tf <- sample_acts_tfs %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>% t()
w <- apply(mat_tf, 1, sd) >= 0.65
mat_tf <- mat_tf[w,]

## ML :: prepare data

neo <- rbind(mat_sig, mat_tf) %>% t()

neo <- neo %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "id") %>% 
  as_tibble() %>% 
  rowwise() %>% 
  mutate(patient = strsplit(id, split = "_")[[1]][1]) %>%
  mutate(y = strsplit(id, split = "_")[[1]][2]) %>%
  select(y, everything()) %>%
  select(-c("id", "patient")) %>% 
  mutate(y = factor(y, levels = c("post", "pre")))

set.seed(123)
neo_split <- initial_split(neo)
neo_train <- training(neo_split)
neo_test <- testing(neo_split)

neo_recipe <- training(neo_split) %>%
  recipe(y ~.) %>%
  step_center(all_predictors(), -all_outcomes()) %>%
  step_scale(all_predictors(), -all_outcomes()) %>%
  prep()

## RF (simple, no tuning)
#https://rviews.rstudio.com/2019/06/19/a-gentle-intro-to-tidymodels/

neo_testing <- neo_recipe %>%
  bake(testing(neo_split)) 

neo_training <- juice(neo_recipe)

neo_ranger <- rand_forest(trees = 100, mode = "classification") %>%
  set_engine("ranger", importance = "permutation") %>%
  fit(y ~ ., data = neo_training)

#predict(neo_ranger, neo_testing)

rf_confusion_mat <- neo_ranger %>%
  predict(neo_testing) %>%
  bind_cols(neo_testing) %>%
  metrics(truth = y, estimate = .pred_class)

outF <- paste(outFolder,"/neo.ML.rf.conf_mat.txt", sep = "")
write_tsv(x = rf_confusion_mat, file = outF)

plt <- neo_ranger %>%
  vip(geom = "col", num_features = 40) + theme_bw()

rf_feature_importance <- vip:::vi(neo_ranger)

outF <- paste(outFolder,"/neo.ML.rf.features.importance.txt", sep = "")
write_tsv(x = rf_feature_importance, file = outF)

#Prepare data for heat map
neo_mat <- rbind(mat_sig, mat_tf)
feat_oi <- rf_feature_importance %>% head(40) %>% pull(Variable)
neo_mat_filt <- neo_mat[feat_oi,]

neo_mat_filt_scaled <- t(apply(neo_mat_filt, 1, scale))
colnames(neo_mat_filt_scaled) <- colnames(neo_mat_filt)

ha <- HeatmapAnnotation(treatment = meta_full[["treatment"]], 
                        response = meta_full[["response"]],
                        col = ha_cols)

HM <- Heatmap(neo_mat_filt_scaled,
              heatmap_legend_param = list(title = "Z-score"),
              col = colorRamp2(c(-2,0,2), c("navy", "white", "red")),
              row_gap = unit(c(3), "mm"),
              row_names_gp = grid::gpar(fontsize = 8),
              column_names_gp = grid::gpar(fontsize = 8),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(3, "cm"),
              border = TRUE,
              top_annotation = ha)

outF <- paste(outFolder,"/neo.ML.rf.features.heatmap.pdf", sep = "")
pdf(outF, width=9, height=7.5)
draw(HM)
dev.off()

## LASSO

lasso_spec <- logistic_reg(penalty = 0.01, mixture = 1) %>%
  set_engine("glmnet")

wf <- workflow() %>%
  add_recipe(neo_recipe)

lasso_fit <- wf %>%
  add_model(lasso_spec) %>%
  fit(data = neo_train)

#feature importance
lasso_feature_importance <- lasso_fit %>%
  pull_workflow_fit() %>%
  tidy()

outF <- paste(outFolder,"/neo.ML.lasso.features.importance.txt", sep = "")
write_tsv(x = lasso_feature_importance, file = outF)

#performance on testing
lasso_confusion_mat <- lasso_fit %>%
  predict(neo_test) %>%
  bind_cols(neo_test) %>%
  mutate(y = as.factor(y)) %>%
  metrics(truth = y, estimate = .pred_class)

outF <- paste(outFolder,"/neo.ML.lasso.conf_mat.txt", sep = "")
write_tsv(x = lasso_confusion_mat, file = outF)

#Prepare data for heat map
neo_mat <- rbind(mat_sig, mat_tf)
feat_oi <- lasso_feature_importance %>% dplyr::filter(estimate != 0 & term != "(Intercept)") %>% pull(term)
neo_mat_filt <- neo_mat[feat_oi,]

neo_mat_filt_scaled <- t(apply(neo_mat_filt, 1, scale))
colnames(neo_mat_filt_scaled) <- colnames(neo_mat_filt)

ha <- HeatmapAnnotation(treatment = meta_full[["treatment"]], 
                        response = meta_full[["response"]],
                        col = ha_cols)

HM <- Heatmap(neo_mat_filt_scaled,
              heatmap_legend_param = list(title = "Z-score"),
              col = colorRamp2(c(-2,0,2), c("navy", "white", "red")),
              row_gap = unit(c(3), "mm"),
              row_names_gp = grid::gpar(fontsize = 8),
              column_names_gp = grid::gpar(fontsize = 8),
              cluster_columns = TRUE,
              cluster_rows = TRUE,
              row_title_rot = 0, 
              column_title_rot = 90,
              row_dend_width = unit(3, "cm"),
              column_dend_height = unit(3, "cm"),
              border = TRUE,
              top_annotation = ha)

outF <- paste(outFolder,"/neo.ML.lasso.features.heatmap.pdf", sep = "")
pdf(outF, width=9, height=4.5)
draw(HM)
dev.off()

