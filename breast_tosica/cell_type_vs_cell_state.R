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

# load comparison pairs
comparisons <- read.xls("../metadata_muw_server.xlsx", sheet = "comparisons") %>%
  filter(untreated %in% c("MCF7_RM", "MCF7-D538G_RM", "MCF7_Lted", "MCF7_TamR"))
  
# calculate fractions
cells_keep <- pred_1 %>%
  filter(Prediction != "Lumsec-prol") %>% pull(cell_id)

fractions_1 <- pred_1 %>%
  filter(Prediction != "Lumsec-prol") %>%    # filter out proliferative fraction
  group_by(cell_type, Prediction) %>%
  summarize(cell_count = n()) %>%
  mutate(fraction = cell_count / sum(cell_count)) %>%
  mutate(Prediction = str_replace_all(Prediction, "-", "_"))

fractions_1_wide <- fractions_1 %>% 
  select(-cell_count) %>% 
  pivot_wider(names_from = Prediction, 
              values_from = fraction, 
              values_fill = 0)

for (line in 1:nrow(comparisons)) {
  untreated <- fractions_1_wide %>% filter(cell_type == comparisons[line,]$untreated)
  treated <- fractions_1_wide %>% filter(cell_type == comparisons[line,]$treated)
  delta <- (treated[1,2:length(treated)] - untreated[1,2:length(untreated)]) %>%
    mutate(cell_type = treated$cell_type)
  
  if (line == 1) {
    deltas <- delta
  } else {
    deltas <- rbind(deltas, delta)
  }
}

deltas_cell_type <- comparisons %>%
  full_join(deltas, by = c("treated" = "cell_type"))

fractions_lumHR_vs_rest <- pred_1 %>% # in this table it's only LumHR vs. the rest
  filter(Prediction != "Lumsec-prol") %>%    # filter out proliferative fraction
  mutate(Prediction_LumHR = ifelse(Prediction %in% c("LummHR-major", "LummHR-active", "LummHR-SCGB"), "LumHR", "other")) %>%
  group_by(cell_type, Prediction_LumHR) %>%
  summarize(cell_count = n()) %>%
  mutate(fraction = cell_count / sum(cell_count)) %>%
  mutate(Prediction_LumHR = str_replace_all(Prediction_LumHR, "-", "_"))

fractions_lumHR_vs_rest_wide <- fractions_lumHR_vs_rest %>% 
  select(-cell_count) %>% 
  pivot_wider(names_from = Prediction_LumHR, 
              values_from = fraction, 
              values_fill = 0)

for (line in 1:nrow(comparisons)) {
  untreated <- fractions_lumHR_vs_rest_wide %>% filter(cell_type == comparisons[line,]$untreated)
  treated <- fractions_lumHR_vs_rest_wide %>% filter(cell_type == comparisons[line,]$treated)
  delta <- (treated[1,2:length(treated)] - untreated[1,2:length(untreated)]) %>%
    mutate(cell_type = treated$cell_type)
  
  if (line == 1) {
    deltas_LumHR <- delta
  } else {
    deltas_LumHR <- rbind(deltas_LumHR, delta)
  }
}

deltas_cell_type_LumHR <- comparisons %>%
  full_join(deltas_LumHR, by = c("treated" = "cell_type")) %>%
  select(-other)
  
# calculate factor usage fractions
NMF_res <- readRDS("../../decipher_full/factor_identification/results/NMF_results_atK.rds")

program_fractions_wide <- NMF_res$MCF7$H %>%
  as.data.frame() %>%
  rownames_to_column(var = "cell") %>%
  filter(cell %in% cells_keep) %>%          # !!! only keep cells that are also used for the prediction deltas
  pivot_longer(cols = starts_with("R33"), names_to = "program", values_to = "value") %>%
  mutate(program = str_replace(program, "R33_Program", "")) %>%
  mutate(cell_type = str_split_i(cell, pattern = "_batch", i = 1),
         larger_zero = ifelse(value > 0, "larger_zero", "zero")) %>%
  group_by(cell_type, program, larger_zero) %>%
  summarize(cell_count = n()) %>%
  mutate(fraction = cell_count / sum(cell_count)) %>%
  filter(larger_zero == "larger_zero") %>%
  select(cell_type, program, fraction) %>%
  pivot_wider(names_from = program, values_from = fraction)

for (line in 1:nrow(comparisons)) {
  untreated <- program_fractions_wide %>% filter(cell_type == comparisons[line,]$untreated)
  treated <- program_fractions_wide %>% filter(cell_type == comparisons[line,]$treated)
  delta <- (treated[1,2:length(treated)] - untreated[1,2:length(untreated)]) %>%
    mutate(cell_type = treated$cell_type)
  
  if (line == 1) {
    deltas_progs <- delta
  } else {
    deltas_progs <- rbind(deltas_progs, delta)
  }
}

deltas_programs <- comparisons %>% 
  full_join(deltas_progs, by = c("treated" = "cell_type"))

# combine cell type and cell state switch information
ct_cs_combined <- deltas_cell_type_LumHR %>% full_join(deltas_programs)

color_table <- read_tsv("../../decipher_full/downstream_analyses/colors_for_umap.txt")
colors <- color_table$color
names(colors) <- color_table$description

plots <- list()
for (program in (deltas_programs %>% select(-treated, -untreated, -group) %>% colnames())) {
  
  first_line <- ct_cs_combined %>% filter(group == "first_line") %>%
    select(treated, group, LumHR, !!program) %>%
    pivot_longer(cols = c("LumHR", program), names_to = "delta_type", values_to = "value") %>%
    mutate(value = ifelse(delta_type == "LumHR", -abs(value), abs(value)))
  
  resistant <- ct_cs_combined %>% filter(group == "resistant") %>%
    select(treated, group, LumHR, !!program) %>%
    pivot_longer(cols = c("LumHR", program), names_to = "delta_type", values_to = "value") %>%
    mutate(value = ifelse(delta_type == "LumHR", -abs(value), abs(value)))
  
  second_line <- ct_cs_combined %>% filter(group == "second_line") %>%
    select(treated, group, LumHR, !!program) %>%
    pivot_longer(cols = c("LumHR", program), names_to = "delta_type", values_to = "value") %>%
    mutate(value = ifelse(delta_type == "LumHR", -abs(value), abs(value)))
  
  first_line_plot <- ggplot(first_line, aes(x = treated, y = value, color = delta_type, fill = treated)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_color_manual(values = c("black", "black")) +
    scale_fill_manual(values = colors) +
    ylim(-0.3, 0.75) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank())
  
  resistant_plot <- ggplot(resistant, aes(x = treated, y = value, color = delta_type, fill = treated)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_color_manual(values = c("black", "black")) +
    scale_fill_manual(values = colors) +
    ylim(-0.3, 0.75) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.title.y=element_blank())
  
  second_line_plot <- ggplot(second_line, aes(x = treated, y = value, color = delta_type, fill = treated)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_color_manual(values = c("black", "black")) +
    scale_fill_manual(values = colors) +
    ylim(-0.3, 0.75) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.y=element_blank()) +
    ylab("Cell type switch fraction (left) vs. cell state switch fraction (right)")
  
  plot_row <- plot_grid(first_line_plot, 
                        resistant_plot,
                        second_line_plot,
                        ncol = 1,
                        align = "v", axis = "bt", rel_heights = c(1, 1, 4.5))
  
  title <- ggdraw() + 
    draw_label(
      paste("Cell type vs. state switch, using gene program ", program, sep = ""),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  
  plots[[program]] <- plot_grid(
    title, plot_row,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  
  pdf(file = paste0("cell_type_vs_state_switch/cell_type_vs_state_lumHRvsRest_method1_program", program, ".pdf"), width = 6, height = 5)
  plot(plots[[program]])
  dev.off()
  
}

#~~~~~~~~~~~~~
# summary plot
##############

ct_cs_filt <- ct_cs_combined %>% 
  filter(!treated %in% c("MCF7_Lted_CDK7i-2d", "MCF7_Lted_Fulv-2d", "MCF7_Lted_Palbo-2d", "MCF7_Lted_Tam-2d", 
                         "MCF7_TamR_CDK7i-2d", "MCF7_TamR_Fulv-2d", "MCF7_TamR_Palbo-2d", "MCF7_TamR_WM-2d")) %>%
  select(treated, untreated, group, LumHR, "14", "23") %>%
  mutate(LumHR = abs(LumHR), `14` = abs(`14`), `23` = abs(`23`))

first_line <- ct_cs_filt %>%
  filter(group == "first_line") %>%
  pivot_longer(cols = c("LumHR", "14", "23"), names_to = "delta_type", values_to = "value")

resistant <- ct_cs_filt %>%
  filter(group == "resistant") %>%
  pivot_longer(cols = c("LumHR", "14", "23"), names_to = "delta_type", values_to = "value")

second_line <- ct_cs_filt %>%
  filter(group == "second_line") %>%
  pivot_longer(cols = c("LumHR", "14", "23"), names_to = "delta_type", values_to = "value")

first_line_plot <- ggplot(first_line, aes(x = treated, y = value, shape = delta_type, fill = treated)) +
  geom_point(size = 4) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c("LumHR" = 21, "14" = 24, "23" = 25)) +
  coord_flip() +
  theme_bw() +
  ylim(0, 0.65) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank())

resistant_plot <- ggplot(resistant, aes(x = treated, y = value, shape = delta_type, fill = treated)) +
  geom_point(size = 4) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c("LumHR" = 21, "14" = 24, "23" = 25)) +
  coord_flip() +
  theme_bw() +
  ylim(0, 0.65) +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank())

second_line_plot <- ggplot(second_line, aes(x = treated, y = value, shape = delta_type, fill = treated)) +
  geom_point(size = 4) +
  scale_fill_manual(values = colors) +
  scale_shape_manual(values = c("LumHR" = 21, "14" = 24, "23" = 25)) +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.y=element_blank()) +
  ylim(0, 0.65) +
  ylab("Cell type vs. state switch (delta fraction of cells)")

summary_plot <- plot_grid(first_line_plot, 
          resistant_plot,
          second_line_plot,
          ncol = 1,
          align = "v", axis = "bt", rel_heights = c(1, 1, 3))

write_tsv(ct_cs_filt, file = "cell_type_vs_state_switch/cell_type_vs_state_lumHRvsRest_method1_summary.tsv")

pdf(file = paste0("cell_type_vs_state_switch/cell_type_vs_state_lumHRvsRest_method1_summary", ".pdf"), width = 5, height = 4)
plot(summary_plot)
dev.off()
  