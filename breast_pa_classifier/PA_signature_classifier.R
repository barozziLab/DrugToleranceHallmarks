## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: PA_signature_classifier.R
##
## Description: 
#~
## Logistic regression classifier for breast data using PA up and down signatures
##
## Authors: 
#~
## Stephan Gruener & Iros Barozzi
##
## License: 
#~
## GNU GPL v3
## Copyright 2024 
## Copyright Iros Barozzi Stephan Gruener
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Notes:
#~
##
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~
# Data preparation
##################

# color scheme
colors <- read.xls("../data/colors_updated_2024.xlsx") %>%
  mutate(description = str_split_i(string = description, pattern = "MCF7.", i = 2))
col <- colors$color_dark
names(col) <- colors$description

# seurat object metadata
data <- read_tsv(file = "../breast_seurat/MCF7_cca_kit_5k_feat/so_metadata_post_integration.tsv")

# metadata
metadata <- read.xls("../data/metadata_filtered_ordered.xlsx") %>%
  mutate(orig.ident = paste0(description, "_batch-", batch)) %>%
  mutate(description = str_split_i(string = description, pattern = "MCF7.", i = 2))

data_noRibo <- data %>% left_join(metadata, by = c("orig.ident" = "orig.ident")) %>%
  select(rowname, description, class, post_integration_pa_up_noRibo1, post_integration_pa_down_noRibo1) %>%
  mutate(description = factor(description, levels = unique(metadata$description)),
         class = factor(class, levels = unique(metadata$class)))

data_ribo <- data %>% left_join(metadata, by = c("orig.ident" = "orig.ident")) %>%
  select(rowname, description, class, post_integration_pa_up1, post_integration_pa_down1) %>%
  mutate(description = factor(description, levels = unique(metadata$description)),
         class = factor(class, levels = unique(metadata$class)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Train classifyer (WM vs. RM)
##############################

full_data <- data_noRibo # select which dataset to train the model with

# Train model only on WM vs. RM
training_full_WM_RM <- full_data %>% 
  filter(description %in% c("RM", "WM-2d"))

pa_up_distr <- ggplot(training_full_WM_RM,
       aes(x = description, y = post_integration_pa_up_noRibo1, fill = description)) +
  geom_violin(position = position_dodge(0.8), scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.05, outlier.shape = NA) +
  scale_fill_brewer(palette = "Accent") +
  theme_classic()

pa_down_distr <- ggplot(training_full_WM_RM,
                      aes(x = description, y = post_integration_pa_down_noRibo1, fill = description)) +
  geom_violin(position = position_dodge(0.8), scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.05, outlier.shape = NA) +
  scale_fill_brewer(palette = "Accent") +
  theme_classic()

pdf("results/pa_up_distr.pdf", width = 4, height = 4)
plot(pa_up_distr)
dev.off()

pdf("results/pa_down_distr.pdf", width = 4, height = 4)
plot(pa_down_distr)
dev.off()

# Reorder the conditions (make factor)
full_data <- full_data %>% mutate(description = factor(description, levels = colors$description))

# Add barplots to visualize the number of observations in each group (also in above plots)
pa_up_all <- ggplot(full_data,
                    aes(x = description, y = post_integration_pa_up_noRibo1, fill = description)) +
  geom_violin(position = position_dodge(0.8), scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.05, outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = col) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))

pa_down_all <- ggplot(full_data,
                      aes(x = description, y = post_integration_pa_down_noRibo1, fill = description)) +
  geom_violin(position = position_dodge(0.8), scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.05, outlier.shape = NA) +
  theme_classic() +
  scale_fill_manual(values = col) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0))

pdf("results/pa_up_all_distr.pdf", width = 7, height = 4)
plot(pa_up_all)
dev.off()

pdf("results/pa_down_all_distr.pdf", width = 7, height = 4)
plot(pa_down_all)
dev.off()

# split training data
set.seed(2503)

training_full_WM_RM <- training_full_WM_RM %>%
  mutate(description = factor(description, levels = c("RM", "WM-2d")))

split_1 <- initial_split(training_full_WM_RM, prop = 0.75, strata = "description")
training_train <- split_1 %>% training()
# match number of observations per class
n_obs <- training_train %>% group_by(description) %>% summarize(count = n()) %>% pull(count) %>% min()
RM_sampled <- training_train %>% filter(description == "RM") %>%
  sample_n(size = n_obs)
training_train_balanced <- training_train %>% filter(description == "WM-2d") %>% rbind(RM_sampled)

# RM and WM ids that can be used for the final analysis
RM_WM_outcome <- split_1 %>% testing()

split_2 <- initial_split(training_train_balanced, prop = 0.75, strata = "description")
training_train_model <- split_2 %>% training()
training_model_eval <- split_2 %>% testing()

# model specification
logistic_model <- logistic_reg(
  mode = "classification",
  engine = "glm")

# fit model
logistic_fit_WM_RM <- fit(logistic_model, description ~ post_integration_pa_up_noRibo1 + post_integration_pa_down_noRibo1,
                    data = training_train_model)

# evaluate model
logistic_class_preds <- predict(logistic_fit_WM_RM, new_data = training_model_eval, type = 'class')
logistic_prob_preds <- predict(logistic_fit_WM_RM, new_data = training_model_eval, type = 'prob')
# bind_cols
logistic_results <- training_model_eval %>% 
  bind_cols(logistic_class_preds, logistic_prob_preds)
# evaluation
logistic_summary <- conf_mat(logistic_results, truth = description, estimate = .pred_class) %>%
  summary()
roc_eval_training <- logistic_results %>% roc_curve(truth = description, `.pred_RM`) %>%
  autoplot()
roc_wm_rm <- logistic_results %>% roc_curve(truth = description, `.pred_RM`) %>% mutate(Model = "WM vs. RM")
roc_auc(logistic_results, truth = description, `.pred_RM`)

pdf("results/roc_curve_training_WM_RM.pdf", width = 3, height = 3)
plot(roc_eval_training)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~
# WM vs. RM predictions
#######################
analysis_data <- full_data %>% filter(!description %in% c("RM", "WM-2d")) %>% rbind(RM_WM_outcome)

model_class_predictions_all_samples <- predict(logistic_fit_WM_RM, new_data = analysis_data, type = 'class')
model_prob_predictions_all_samples <- predict(logistic_fit_WM_RM, new_data = analysis_data, type = 'prob')
model_pred_results_all_samples <- analysis_data %>%
  bind_cols(model_class_predictions_all_samples, model_prob_predictions_all_samples) %>%
  mutate(description = factor(description, levels = rev(colors$description)))

stats_wm_rm <- model_pred_results_all_samples %>% group_by(description) %>%
  summarize(n_cells = n(),
            n_gt.5 = sum(`.pred_WM-2d` > 0.5),
            p_gt.5 = n_gt.5 / n_cells,
            avg_score = mean(`.pred_WM-2d`))
  
model_pred_results_all_samples <- model_pred_results_all_samples %>%
  left_join(stats_wm_rm, by = "description")

pred_scores <- ggplot(model_pred_results_all_samples, aes(x = `.pred_WM-2d`, y = description, fill = description)) +
  geom_violin(position = position_dodge(0.8), scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.2, outlier.shape = NA, fill = "white") +
  theme_classic() +
  scale_fill_manual(values = col) +
  theme(legend.position = "none") +
  #theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  xlab("prediction score: WM-2d") +
  ylab("") +
  theme(axis.text.y = element_text(hjust = 0.5)) +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, -16), "pt"))

pred_scores_stat <- ggplot(model_pred_results_all_samples, aes(x = p_gt.5, y = description, fill = description)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  scale_fill_manual(values = col) +
  theme(legend.position = "none") +
  #theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  xlab("percentage of cells > 0.5") +
  scale_x_reverse() +
  scale_y_discrete(position = "right") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

rm_wm_grid <- cowplot::plot_grid(pred_scores_stat, pred_scores, 
                   rel_widths = c(0.6, 1),
                   ncol = 2)

pdf("results/prediction_scores_WM_RM.pdf", width = 5, height = 5)
plot(rm_wm_grid)
dev.off()

stats_wm_rm <- stats_wm_rm %>% dplyr::rename(sample = description)
write_tsv(stats_wm_rm, "results/wm_rm_stats.tsv")

#~~~~~~~~~~~
# WM vs. TAM
############

# Train model only on WM vs. TAM
training_full_WM_TAM <- full_data %>% 
  filter(description %in% c("Tam-2d", "WM-2d"))

# split training data
set.seed(2503)

training_full_WM_TAM <- training_full_WM_TAM %>%
  mutate(description = factor(description, levels = c("Tam-2d", "WM-2d")))

split_1 <- initial_split(training_full_WM_TAM, prop = 0.75, strata = "description")
training_train <- split_1 %>% training()
# match number of observations per class
n_obs <- training_train %>% group_by(description) %>% summarize(count = n()) %>% pull(count) %>% min()
TAM_sampled <- training_train %>% filter(description == "Tam-2d") %>%
  sample_n(size = n_obs)
training_train_balanced <- training_train %>% filter(description == "WM-2d") %>% rbind(TAM_sampled)

# RM and WM ids that can be used for the final analysis
TAM_WM_outcome <- split_1 %>% testing()

split_2 <- initial_split(training_train_balanced, prop = 0.75, strata = "description")
training_train_model <- split_2 %>% training()
training_model_eval <- split_2 %>% testing()

# model specification
logistic_model <- logistic_reg(
  mode = "classification",
  engine = "glm")

# fit model
logistic_fit_WM_TAM <- fit(logistic_model, description ~ post_integration_pa_up_noRibo1 + post_integration_pa_down_noRibo1,
                          data = training_train_model)

# evaluate model
logistic_class_preds <- predict(logistic_fit_WM_TAM, new_data = training_model_eval, type = 'class')
logistic_prob_preds <- predict(logistic_fit_WM_TAM, new_data = training_model_eval, type = 'prob')
# bind_cols
logistic_results <- training_model_eval %>% 
  bind_cols(logistic_class_preds, logistic_prob_preds)
# evaluation
logistic_summary <- conf_mat(logistic_results, truth = description, estimate = .pred_class) %>%
  summary()
roc_eval_training <- logistic_results %>% roc_curve(truth = description, `.pred_Tam-2d`) %>%
  autoplot()
roc_wm_tam <- logistic_results %>% roc_curve(truth = description, `.pred_Tam-2d`) %>% mutate(Model = "WM vs. Tam")
roc_auc(logistic_results, truth = description, `.pred_Tam-2d`)

pdf("results/roc_curve_training_WM_TAM.pdf", width = 4, height = 4)
plot(roc_eval_training)
dev.off()

# combined roc curve
roc_comb <- bind_rows(roc_wm_rm, roc_wm_tam)
roc_col <- col[c("RM", "Tam-2d")]
names(roc_col) <- c("WM vs. RM", "WM vs. Tam")

roc_comb_plt <- 
  ggplot(roc_comb, aes(x = 1 - specificity, y = sensitivity, color = Model)) +
  geom_line(linewidth = 1.5) +
  theme_bw() +
  scale_color_manual(values = roc_col) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7, 0.6),
        legend.background = element_rect(fill = "white",
                                         color = "black",
                                         linewidth = 0.2))

pdf("results/roc_curve_training_combined.pdf", width = 3, height = 2.9)
plot(roc_comb_plt)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# WM vs. TAM model predictions
##############################
analysis_data <- full_data %>% filter(!description %in% c("Tam-2d", "WM-2d")) %>% rbind(TAM_WM_outcome)

model_class_predictions_all_samples <- predict(logistic_fit_WM_TAM, new_data = analysis_data, type = 'class')
model_prob_predictions_all_samples <- predict(logistic_fit_WM_TAM, new_data = analysis_data, type = 'prob')
model_pred_results_all_samples <- analysis_data %>%
  bind_cols(model_class_predictions_all_samples, model_prob_predictions_all_samples)


pred_scores <- ggplot(model_pred_results_all_samples, aes(x = description, y = `.pred_WM-2d`, fill = description)) +
  geom_violin(position = position_dodge(0.8), scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.05, outlier.shape = NA, fill = "white") +
  theme_classic() +
  scale_fill_manual(values = col) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
  geom_hline(yintercept = 0.5)

pdf("results/prediction_scores_WM_TAM.pdf", width = 7, height = 4)
plot(pred_scores)
dev.off()

