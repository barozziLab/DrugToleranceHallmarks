# load data and metadata (and join)
metadata <- read.xls("../metadata_filtered_ordered.xlsx", sheet = "filtered") %>%
  as_tibble() %>%
  mutate(orig.ident = paste0(description, "_batch-", batch))

pred_1 <- read_tsv("predictions_full_1.txt") %>%
  select(-cell_type_id) %>%
  left_join(metadata %>% select(-c("kit", "batch")), by = c("orig.ident" = "orig.ident"))
colnames(pred_1)[1] <- "cell_id"
pred_2 <- read_tsv("predictions_full_2.txt") %>%
  select(-cell_type_id) %>%
  left_join(metadata %>% select(-c("kit", "batch")), by = c("orig.ident" = "orig.ident"))
colnames(pred_2)[1] <- "cell_id"
pred_3 <- read_tsv("predictions_full_3.txt") %>%
  left_join(pred_2 %>% select(-c("Prediction", "Probability")), by = c("...1" = "cell_id")) # since the orig.ident column was not saved for this method, join with metadata from method 2
colnames(pred_3)[1] <- "cell_id"

# load cell-cycle annotation (and join)
cc_anno <- read_tsv("../metadata_MCF7_full_cc.tsv") %>% select(cell_id, Phase)

pred_1 <- pred_1 %>%
  left_join(cc_anno, by = c("cell_id" = "cell_id"))
pred_2 <- pred_2 %>%
  left_join(cc_anno, by = c("cell_id" = "cell_id"))
pred_3 <- pred_3 %>%
  left_join(cc_anno, by = c("cell_id" = "cell_id"))

# order of conditions
order <- metadata$description %>% unique()

# define prediction colors
prediction_colors <- c("basal" = "#970001", 
                       "LummHR-major" = "#FFBDDA", 
                       "LummHR-active" = "#FF515A", 
                       "LummHR-SCGB" = "#FE01FF",
                       "Lumsec-major" = "#A16664",
                       "Lumsec-basal" = "#FFD024",
                       "Lumsec-HLA" = "#FF8000",
                       "Lumsec-KIT" = "#BC8A00",
                       "Lumsec-lac" = "#ECCCB8",
                       "Lumsec-myo" = "#B7B900",
                       "Lumsec-prol" = "#3D4E22")

prediction_colors_broad <- c("Basal" = "#970001", 
                       "LumHR" = "#FFBDDA",
                       "LumSec" = "#3D4E22")

# predictions by treatment
p1_numbers <- pred_1 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 1") +
  theme(axis.title.y = element_blank())

p2_numbers <- pred_2 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 2") +
  theme(axis.title.y = element_blank())

p3_numbers <- pred_3 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 3") +
  theme(axis.title.y = element_blank())

cc_numbers <- pred_1 %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  #scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Cell Cycle") +
  theme(axis.title.y = element_blank())

p1_perc <- pred_1 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 1") +
  theme(axis.title.y = element_blank())

p2_perc <- pred_2 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 2") +
  theme(axis.title.y = element_blank())

p3_perc <- pred_3 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 3") +
  theme(axis.title.y = element_blank())

cc_perc <- pred_1 %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Phase)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  #scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Cell Cycle") +
  theme(axis.title.y = element_blank())

fractions_ct_1 <- pred_1 %>%
  mutate(Prediction = factor(Prediction)) %>%
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>% 
  mutate('Percentage_Lumsec-prol' = count / sum(count)) %>%
  filter(Prediction == "Lumsec-prol")

fractions_cc_1 <- pred_1 %>%
  mutate(Phase = factor(Phase)) %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>% 
  mutate(Percentage_G1 = count / sum(count)) %>%
  filter(Phase == "G1")

fractions_1 <- full_join(fractions_ct_1, fractions_cc_1, by = c("cell_type" = "cell_type")) %>% mutate(method = "method1")
frac_plt_1 <- ggplot(fractions_1, aes(x = Percentage_G1, y = `Percentage_Lumsec-prol`, color = cell_type)) +
  geom_point() +
  geom_text_repel(aes(label = cell_type), size = 3) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Prediction method 1") +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", aes(x = Percentage_G1, y = `Percentage_Lumsec-prol`), inherit.aes = F) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.5, cor.coef.name = "rho", aes(x = Percentage_G1, y = `Percentage_Lumsec-prol`), inherit.aes = F)

fractions_ct_2 <- pred_2 %>%
  mutate(Prediction = factor(Prediction)) %>%
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>% 
  mutate('Percentage_Lumsec-prol' = count / sum(count)) %>%
  filter(Prediction == "Lumsec-prol")

fractions_cc_2 <- pred_2 %>%
  mutate(Phase = factor(Phase)) %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>% 
  mutate(Percentage_G1 = count / sum(count)) %>%
  filter(Phase == "G1")

fractions_2 <- full_join(fractions_ct_2, fractions_cc_2, by = c("cell_type" = "cell_type")) %>% mutate(method = "method2")
frac_plt_2 <- ggplot(fractions_2, aes(x = Percentage_G1, y = `Percentage_Lumsec-prol`, col = cell_type)) +
  geom_point() +
  geom_text_repel(aes(label = cell_type), size = 3) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Prediction method 2") +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", aes(x = Percentage_G1, y = `Percentage_Lumsec-prol`), inherit.aes = F) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.5, cor.coef.name = "rho", aes(x = Percentage_G1, y = `Percentage_Lumsec-prol`), inherit.aes = F)

fractions_ct_3 <- pred_3 %>%
  mutate(Prediction = factor(Prediction)) %>%
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>% 
  mutate('Percentage_Lumsec-prol' = count / sum(count)) %>%
  filter(Prediction == "Lumsec-prol")

fractions_cc_3 <- pred_3 %>%
  mutate(Phase = factor(Phase)) %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>% 
  mutate(Percentage_G1 = count / sum(count)) %>%
  filter(Phase == "G1")

fractions_3 <- full_join(fractions_ct_3, fractions_cc_3, by = c("cell_type" = "cell_type")) %>% mutate(method = "method3")
frac_plt_3 <- ggplot(fractions_3, aes(x = Percentage_G1, y = `Percentage_Lumsec-prol`, col = cell_type)) +
  geom_point() +
  geom_text_repel(aes(label = cell_type), size = 3) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Prediction method 3") +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", aes(x = Percentage_G1, y = `Percentage_Lumsec-prol`), inherit.aes = F) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.5, cor.coef.name = "rho", aes(x = Percentage_G1, y = `Percentage_Lumsec-prol`), inherit.aes = F)

pred_grid <- plot_grid(p1_numbers, p2_numbers, p3_numbers, cc_numbers, p1_perc, p2_perc, p3_perc, cc_perc, ncol = 4)
frac_grid <- plot_grid(frac_plt_1, frac_plt_2, frac_plt_3, ncol = 3)

fractions_summary <- rbind(fractions_1, fractions_2, fractions_3)
write_tsv(fractions_summary, "fractions_prol_cc_summary.tsv")

predictions_summary_data <- rbind(
  p1_numbers$data %>% mutate(method = "method1"),
  p2_numbers$data %>% mutate(method = "method2"),
  p3_numbers$data %>% mutate(method = "method3")
) %>% select(-.group)
predictions_summary_cc_data <- cc_numbers$data %>% select(-.group)
write_tsv(predictions_summary_data, "predictions_summary_data.tsv")
write_tsv(predictions_summary_cc_data, "cc_data_full.tsv")

ggsave(filename = "predictions_summary.pdf", plot = pred_grid, device = "pdf", width = 22, height = 10)
ggsave(filename = "predictions_fractions_summary.pdf", plot = frac_grid, device = "pdf", width = 19, height = 7)

#~~~~~~~~~~~~~~~~~~~~~~~~~
# broad cell type analysis
##########################

pred_1_b <- read_tsv("predictions_hr_sec_1.txt") %>%
  select(-cell_type_id) %>%
  left_join(metadata %>% select(-c("kit", "batch")), by = c("orig.ident" = "orig.ident"))
colnames(pred_1_b)[1] <- "cell_id"
pred_2_b <- read_tsv("predictions_hr_sec_2.txt") %>%
  select(-cell_type_id) %>%
  left_join(metadata %>% select(-c("kit", "batch")), by = c("orig.ident" = "orig.ident"))
colnames(pred_2_b)[1] <- "cell_id"
pred_3_b <- read_tsv("predictions_hr_sec_3.txt") %>%
  left_join(pred_2_b %>% select(-c("Prediction", "Probability")), by = c("...1" = "cell_id")) # since the orig.ident column was not saved for this method, join with metadata from method 2
colnames(pred_3_b)[1] <- "cell_id"

pred_1_b <- pred_1_b %>%
  left_join(cc_anno, by = c("cell_id" = "cell_id"))
pred_2_b <- pred_2_b %>%
  left_join(cc_anno, by = c("cell_id" = "cell_id"))
pred_3_b <- pred_3_b %>%
  left_join(cc_anno, by = c("cell_id" = "cell_id"))

# predictions by treatment
p1_b_numbers <- pred_1_b %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors_broad) +
  coord_flip() +
  ggtitle("Predictions method 1") +
  theme(axis.title.y = element_blank())

p2_b_numbers <- pred_2_b %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors_broad) +
  coord_flip() +
  ggtitle("Predictions method 2") +
  theme(axis.title.y = element_blank())

p3_b_numbers <- pred_3_b %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors_broad) +
  coord_flip() +
  ggtitle("Predictions method 3") +
  theme(axis.title.y = element_blank())

p1_b_perc <- pred_1_b %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors_broad) +
  coord_flip() +
  ggtitle("Predictions method 1") +
  theme(axis.title.y = element_blank())

p2_b_perc <- pred_2_b %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors_broad) +
  coord_flip() +
  ggtitle("Predictions method 2") +
  theme(axis.title.y = element_blank())

p3_b_perc <- pred_3_b %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors_broad) +
  coord_flip() +
  ggtitle("Predictions method 3") +
  theme(axis.title.y = element_blank())

fractions_ct_1_b <- pred_1_b %>%
  mutate(Prediction = factor(Prediction)) %>%
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>% 
  mutate('Percentage_LumSec' = count / sum(count)) %>%
  filter(Prediction == "LumSec")

fractions_cc_1_b <- pred_1_b %>%
  mutate(Phase = factor(Phase)) %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>% 
  mutate(Percentage_G1 = count / sum(count)) %>%
  filter(Phase == "G1")

fractions_1_b <- full_join(fractions_ct_1_b, fractions_cc_1_b, by = c("cell_type" = "cell_type"))
frac_plt_1_b <- ggplot(fractions_1_b, aes(x = Percentage_G1, y = `Percentage_LumSec`, col = cell_type)) +
  geom_point() +
  geom_text_repel(aes(label = cell_type), size = 3) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Prediction method 1") +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", aes(x = Percentage_G1, y = `Percentage_LumSec`), inherit.aes = F) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.5, cor.coef.name = "rho", aes(x = Percentage_G1, y = `Percentage_LumSec`), inherit.aes = F)

fractions_ct_2_b <- pred_2_b %>%
  mutate(Prediction = factor(Prediction)) %>%
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>% 
  mutate('Percentage_LumSec' = count / sum(count)) %>%
  filter(Prediction == "LumSec")

fractions_cc_2_b <- pred_2_b %>%
  mutate(Phase = factor(Phase)) %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>% 
  mutate(Percentage_G1 = count / sum(count)) %>%
  filter(Phase == "G1")

fractions_2_b <- full_join(fractions_ct_2_b, fractions_cc_2_b, by = c("cell_type" = "cell_type"))
frac_plt_2_b <- ggplot(fractions_2_b, aes(x = Percentage_G1, y = `Percentage_LumSec`, col = cell_type)) +
  geom_point() +
  geom_text_repel(aes(label = cell_type), size = 3) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Prediction method 2") +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", aes(x = Percentage_G1, y = `Percentage_LumSec`), inherit.aes = F) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.5, cor.coef.name = "rho", aes(x = Percentage_G1, y = `Percentage_LumSec`), inherit.aes = F)

fractions_ct_3_b <- pred_3_b %>%
  mutate(Prediction = factor(Prediction)) %>%
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>% 
  mutate('Percentage_LumSec' = count / sum(count)) %>%
  filter(Prediction == "LumSec")

fractions_cc_3_b <- pred_3_b %>%
  mutate(Phase = factor(Phase)) %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>% 
  mutate(Percentage_G1 = count / sum(count)) %>%
  filter(Phase == "G1")

fractions_3_b <- full_join(fractions_ct_3_b, fractions_cc_3_b, by = c("cell_type" = "cell_type"))
frac_plt_3_b <- ggplot(fractions_3_b, aes(x = Percentage_G1, y = `Percentage_LumSec`, col = cell_type)) +
  geom_point() +
  geom_text_repel(aes(label = cell_type), size = 3) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("Prediction method 3") +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", aes(x = Percentage_G1, y = `Percentage_LumSec`), inherit.aes = F) +
  ggpubr::stat_cor(method = "spearman", label.x = 0.5, cor.coef.name = "rho", aes(x = Percentage_G1, y = `Percentage_LumSec`), inherit.aes = F)

pred_grid_b <- plot_grid(p1_b_numbers, p2_b_numbers, p3_b_numbers, cc_numbers, p1_b_perc, p2_b_perc, p3_b_perc, cc_perc, ncol = 4)
frac_grid_b <- plot_grid(frac_plt_1_b, frac_plt_2_b, frac_plt_3_b, ncol = 3)

predictions_b_summary_data <- rbind(
  p1_b_numbers$data %>% mutate(method = "method1"),
  p2_b_numbers$data %>% mutate(method = "method2"),
  p3_b_numbers$data %>% mutate(method = "method3")
) %>% select(-.group)

write_tsv(predictions_b_summary_data, "predictions_broad_summary_data.tsv")

ggsave(filename = "predictions_broad_summary.pdf", plot = pred_grid_b, device = "pdf", width = 20, height = 10)
ggsave(filename = "predictions_fractions_broad_summary.pdf", plot = frac_grid_b, device = "pdf", width = 19, height = 7)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# same analysis after removing Lumsec-prol
##########################################
pred_1 <- pred_1 %>% filter(Prediction != "Lumsec-prol")
pred_2 <- pred_2 %>% filter(Prediction != "Lumsec-prol")
pred_3 <- pred_3 %>% filter(Prediction != "Lumsec-prol")
pred_1_b <- pred_1_b %>% filter(Prediction != "LumSec")
pred_2_b <- pred_2_b %>% filter(Prediction != "LumSec")
pred_3_b <- pred_3_b %>% filter(Prediction != "LumSec")

# predictions by treatment
p1_numbers_filt <- pred_1 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 1") +
  theme(axis.title.y = element_blank())

p2_numbers_filt <- pred_2 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 2") +
  theme(axis.title.y = element_blank())

p3_numbers_filt <- pred_3 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 3") +
  theme(axis.title.y = element_blank())

cc_numbers_filt <- pred_1 %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Phase)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  #scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Cell Cycle") +
  theme(axis.title.y = element_blank())

p1_perc_filt <- pred_1 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 1") +
  theme(axis.title.y = element_blank())

# p1_perc_filt_poster <- pred_1 %>% 
#   group_by(cell_type, Prediction) %>% 
#   summarize(count = n()) %>%
#   filter(cell_type %in% c("MCF7_RM", "MCF7_WM-2d", "MCF7_Tam-2d", "MCF7_Fulv-2d", "MCF7_TamR", "MCF7_Lted", "MCF7_FulvR",
#                           "MCF7_TamR_WM-7d", "MCF7_TamR_Palbo-7d", "MCF7_TamR_Fulv-7d", "MCF7_TamR_CDK7i-7d",
#                           "MCF7_Lted_Tam-7d", "MCF7_Lted_Palbo-7d", "MCF7_Lted_Fulv-7d", "MCF7_Lted_CDK7i-7d")) %>%
#   mutate(cell_type = factor(cell_type, levels = rev(c("MCF7_RM", "MCF7_WM-2d", "MCF7_Tam-2d", "MCF7_Fulv-2d", "MCF7_TamR", "MCF7_Lted", "MCF7_FulvR",
#                                                   "MCF7_TamR_WM-7d", "MCF7_TamR_Palbo-7d", "MCF7_TamR_Fulv-7d", "MCF7_TamR_CDK7i-7d",
#                                                   "MCF7_Lted_Tam-7d", "MCF7_Lted_Palbo-7d", "MCF7_Lted_Fulv-7d", "MCF7_Lted_CDK7i-7d")))) %>%
#   ggplot(aes(x = cell_type, y = count, fill = Prediction)) +
#   geom_bar(color = "black", stat = "identity", position = "fill") +
#   theme_bw() +
#   scale_fill_manual(values = prediction_colors) +
#   coord_flip() +
#   theme(axis.title.y = element_blank())

p2_perc_filt <- pred_2 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 2") +
  theme(axis.title.y = element_blank())

p3_perc_filt <- pred_3 %>% 
  group_by(cell_type, Prediction) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Predictions method 3") +
  theme(axis.title.y = element_blank())

cc_perc_filt <- pred_1 %>%
  group_by(cell_type, Phase) %>% 
  summarize(count = n()) %>%
  ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Phase)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  #scale_fill_manual(values = prediction_colors) +
  coord_flip() +
  ggtitle("Cell Cycle") +
  theme(axis.title.y = element_blank())

pred_grid_filt <- plot_grid(p1_numbers_filt, p2_numbers_filt, p3_numbers_filt, cc_numbers_filt, p1_perc_filt, p2_perc_filt, p3_perc_filt, cc_perc_filt, ncol = 4)

predictions_summary_data_filt <- rbind(
  p1_numbers_filt$data %>% mutate(method = "method1"),
  p2_numbers_filt$data %>% mutate(method = "method2"),
  p3_numbers_filt$data %>% mutate(method = "method3")
) %>% select(-.group)
predictions_summary_cc_data_filt <- cc_numbers_filt$data %>% select(-.group)
write_tsv(predictions_summary_data_filt, "predictions_summary_data_filt_no_prol.tsv")
write_tsv(predictions_summary_cc_data_filt, "cc_data_filt_no_prol.tsv")

ggsave(filename = "predictions_summary_no_lumsec_prol.pdf", plot = pred_grid_filt, device = "pdf", width = 22, height = 10)
ggsave(filename = "p1_perc_filt_poster.pdf", plot = p1_perc_filt_poster, device = "pdf", width = 5, height = 4)

# # predictions by treatment (broad annotation)
# p1_b_numbers_filt <- pred_1_b %>% 
#   group_by(cell_type, Prediction) %>% 
#   summarize(count = n()) %>%
#   ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
#   geom_bar(stat = "identity", position = "stack") +
#   theme_bw() +
#   scale_fill_manual(values = prediction_colors_broad) +
#   coord_flip() +
#   ggtitle("Predictions method 1") +
#   theme(axis.title.y = element_blank())
# 
# p2_b_numbers_filt <- pred_2_b %>% 
#   group_by(cell_type, Prediction) %>% 
#   summarize(count = n()) %>%
#   ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
#   geom_bar(stat = "identity", position = "stack") +
#   theme_bw() +
#   scale_fill_manual(values = prediction_colors_broad) +
#   coord_flip() +
#   ggtitle("Predictions method 2") +
#   theme(axis.title.y = element_blank())
# 
# p3_b_numbers_filt <- pred_3_b %>% 
#   group_by(cell_type, Prediction) %>% 
#   summarize(count = n()) %>%
#   ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
#   geom_bar(stat = "identity", position = "stack") +
#   theme_bw() +
#   scale_fill_manual(values = prediction_colors_broad) +
#   coord_flip() +
#   ggtitle("Predictions method 3") +
#   theme(axis.title.y = element_blank())
# 
# p1_b_perc_filt <- pred_1_b %>% 
#   group_by(cell_type, Prediction) %>% 
#   summarize(count = n()) %>%
#   ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
#   geom_bar(stat = "identity", position = "fill") +
#   theme_bw() +
#   scale_fill_manual(values = prediction_colors_broad) +
#   coord_flip() +
#   ggtitle("Predictions method 1") +
#   theme(axis.title.y = element_blank())
# 
# p2_b_perc_filt <- pred_2_b %>% 
#   group_by(cell_type, Prediction) %>% 
#   summarize(count = n()) %>%
#   ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
#   geom_bar(stat = "identity", position = "fill") +
#   theme_bw() +
#   scale_fill_manual(values = prediction_colors_broad) +
#   coord_flip() +
#   ggtitle("Predictions method 2") +
#   theme(axis.title.y = element_blank())
# 
# p3_b_perc_filt <- pred_3_b %>% 
#   group_by(cell_type, Prediction) %>% 
#   summarize(count = n()) %>%
#   ggplot(aes(x = factor(cell_type, levels = rev(order)), y = count, fill = Prediction)) +
#   geom_bar(stat = "identity", position = "fill") +
#   theme_bw() +
#   scale_fill_manual(values = prediction_colors_broad) +
#   coord_flip() +
#   ggtitle("Predictions method 3") +
#   theme(axis.title.y = element_blank())
# 
# pred_grid_b_filt <- plot_grid(p1_b_numbers_filt, p2_b_numbers_filt, p3_b_numbers_filt, cc_numbers_filt, p1_b_perc_filt, p2_b_perc_filt, p3_b_perc_filt, cc_perc_filt, ncol = 4)
# 
# ggsave(filename = "predictions_broad_summary_no_lumsec_prol.pdf", plot = pred_grid_b_filt, device = "pdf", width = 20, height = 10)