## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: run.R
##
## Description: 
#~
## Investigate PA ~ cell cycle correlation
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


library(tidyverse)
library(gdata)
library(Seurat)
library(msigdbr)
library(ComplexHeatmap)

out_folder <- "cc_results/"
cmd <- paste0("mkdir ", out_folder)
system(cmd)


#~~~~~~~~~~~~~~~~~~
#Metadata (Samples)
###################

metadata_ordered <- read.xls("../data/metadata_filtered_ordered.xlsx") %>% 
  as_tibble() %>%
  mutate(sample = paste0(description, "_batch-", batch)) %>%
  mutate(class = ifelse(ESR1_status == "mut", "ESR1_mutant", class)) %>%
  mutate(class = factor(class, levels = c("treatment-naive", "ESR1_mutant", "first_line", "resistant", "second_line"))) %>%
  dplyr::rename(cell_type = description)

#add latest colors
new_cols <- read_csv("../data/colors_for_umap.txt")
new_cols <- new_cols %>% 
  dplyr::rename(cell_type = description) %>% 
  dplyr::rename(treatment_2 = treatment)
metadata_ordered <- metadata_ordered %>%
  left_join(new_cols, by = "cell_type")

coi <- read_tsv("coi.short.txt", col_names = FALSE) %>% pull(1)
coi_all <- read_tsv("coi.long.txt", col_names = FALSE) %>% pull(1)

#prepare color table for coi_all
cols_table <- metadata_ordered %>% 
  select(cell_type, color) %>% 
  unique() %>% 
  column_to_rownames("cell_type")
cols_table <- cols_table[rev(coi_all),,drop=FALSE]


#~~~~~~~~~~~~~~~~
#Data Preparation
#################

#metadata
meta <- read_tsv("so_metadata_post_integration.tsv.gz")

#fix colnames (strip 1s)
colnames(meta) <- gsub("1$", "", colnames(meta))

#filter for coi_all
meta <- meta %>% dplyr::filter(cell_type %in% coi_all)

#define groups of conditions
groups <- read_tsv("groups.txt")
meta <- meta %>% left_join(groups, by = "cell_type")

#re-order factor levels of cell_type
meta <- meta %>% mutate(cell_type = factor(cell_type, levels = rev(coi_all)))

#strip MCF7
meta <- meta %>% mutate(cell_type_2 = gsub("MCF7_|MCF7-", "", cell_type))
meta <- meta %>% mutate(cell_type_2 = factor(cell_type_2, levels = gsub("MCF7_|MCF7-", "", rev(coi_all))))


#~~~~~~~
#PA ~ CC
########

##Hallmarks and PA use post integration --> see plots

##Ribo Yes

qc_plts <- list()

#pre

cor_res <- cor.test(meta$prior_integration_pa_up, meta$prior_integration_pa_down, method = "spearman")
plt_main <- paste0("SCC =", formatC(cor_res$estimate, 3))

qc_plts[["prior_scatter"]] <- ggplot(data = meta, mapping = aes(x = prior_integration_pa_down, y = prior_integration_pa_up, col = Phase)) +
  geom_point(shape = 1, size = 1) +
  ggtitle(plt_main) +
  theme_bw()

qc_plts[["prior_scatter_facet"]] <- ggplot(data = meta, mapping = aes(x = prior_integration_pa_down, y = prior_integration_pa_up, col = Phase)) +
  geom_point(shape = 1, size = 1) +
  facet_wrap(~orig.ident) +
  theme_bw()

out_path <- paste0(out_folder, "/qc_prior_scatter.ribo_yes.png")
png(out_path, width = 600, height = 500, res = 150)
plot(qc_plts[["prior_scatter"]])
dev.off()

out_path <- paste0(out_folder, "/qc_prior_scatter_facet.ribo_yes.png")
png(out_path, width = 1800, height = 1800, res = 150)
plot(qc_plts[["prior_scatter_facet"]])
dev.off()

#post

cor_res <- cor.test(meta$post_integration_pa_up, meta$post_integration_pa_down, method = "spearman")
plt_main <- paste0("SCC =", formatC(cor_res$estimate, 3))

qc_plts[["post_scatter"]] <- ggplot(data = meta, mapping = aes(x = post_integration_pa_down, y = post_integration_pa_up, col = Phase)) +
  geom_point(shape = 1, size = 1) +
  ggtitle(plt_main) +
  theme_bw()

qc_plts[["post_scatter_facet"]] <- ggplot(data = meta, mapping = aes(x = post_integration_pa_down, y = post_integration_pa_up, col = Phase)) +
  geom_point(shape = 1, size = 1) +
  facet_wrap(~orig.ident) +
  theme_bw()

out_path <- paste0(out_folder, "/qc_post_scatter.ribo_yes.png")
png(out_path, width = 600, height = 500, res = 150)
plot(qc_plts[["post_scatter"]])
dev.off()

out_path <- paste0(out_folder, "/qc_post_scatter_facet.ribo_yes.png")
png(out_path, width = 1800, height = 1800, res = 150)
plot(qc_plts[["post_scatter_facet"]])
dev.off()

##Ribo No

qc_plts <- list()

#pre

cor_res <- cor.test(meta$prior_integration_pa_up_noRibo, meta$prior_integration_pa_down_noRibo, method = "spearman")
plt_main <- paste0("SCC =", formatC(cor_res$estimate, 3))

qc_plts[["prior_scatter"]] <- ggplot(data = meta, mapping = aes(x = prior_integration_pa_down_noRibo, y = prior_integration_pa_up_noRibo, col = Phase)) +
  geom_point(shape = 1, size = 1) +
  ggtitle(plt_main) +
  theme_bw()

qc_plts[["prior_scatter_facet"]] <- ggplot(data = meta, mapping = aes(x = prior_integration_pa_down_noRibo, y = prior_integration_pa_up_noRibo, col = Phase)) +
  geom_point(shape = 1, size = 1) +
  facet_wrap(~orig.ident) +
  theme_bw()

out_path <- paste0(out_folder, "/qc_prior_scatter.ribo_no.png")
png(out_path, width = 600, height = 500, res = 150)
plot(qc_plts[["prior_scatter"]])
dev.off()

out_path <- paste0(out_folder, "/qc_prior_scatter_facet.ribo_no.png")
png(out_path, width = 1800, height = 1800, res = 150)
plot(qc_plts[["prior_scatter_facet"]])
dev.off()

#post

cor_res <- cor.test(meta$post_integration_pa_up_noRibo, meta$post_integration_pa_down_noRibo, method = "spearman")
plt_main <- paste0("SCC =", formatC(cor_res$estimate, 3))

qc_plts[["post_scatter"]] <- ggplot(data = meta, mapping = aes(x = post_integration_pa_down_noRibo, y = post_integration_pa_up_noRibo, col = Phase)) +
  geom_point(shape = 1, size = 1) +
  ggtitle(plt_main) +
  theme_bw()

qc_plts[["post_scatter_facet"]] <- ggplot(data = meta, mapping = aes(x = post_integration_pa_down_noRibo, y = post_integration_pa_up_noRibo, col = Phase)) +
  geom_point(shape = 1, size = 1) +
  facet_wrap(~orig.ident) +
  theme_bw()

out_path <- paste0(out_folder, "/qc_post_scatter.ribo_no.png")
png(out_path, width = 600, height = 500, res = 150)
plot(qc_plts[["post_scatter"]])
dev.off()

out_path <- paste0(out_folder, "/qc_post_scatter_facet.ribo_no.png")
png(out_path, width = 1800, height = 1800, res = 150)
plot(qc_plts[["post_scatter_facet"]])
dev.off()


#~~~~~~~
#PA ~ CC
########

#PA, post
w <- intersect(grep("post_", colnames(meta)), grep("pa_", colnames(meta)))
sig_names <- colnames(meta)[w]

#Median per PA signature, CC phase
meta_stats <- meta %>% 
  select(cell_type_2, Phase, sig_names[1:length(sig_names)]) %>% 
  pivot_longer(cols = -c("cell_type_2", "Phase")) %>% 
  group_by(cell_type_2, Phase, name) %>% 
  summarise(median = median(value)) %>% 
  arrange(name) %>%
  ungroup()

cc_plts <- c()

for (sig_name in sig_names) {

  cc_plts[[sig_name]] <- ggplot(meta, aes(x=cell_type_2, y=get(sig_name), fill=cell_type_2)) +
    geom_boxplot() + 
    scale_fill_manual(values = cols_table[,1]) + 
    facet_wrap(~Phase) +
    coord_flip() +
    ylim(-0.75, 0.75) +
    ylab(sig_name) +
    #geom_hline(yintercept = 0, color = "black") +
    theme_bw() + 
    theme(legend.position = "none")

}

out_path <- paste0(out_folder, "/PA-scores.boxplots.pdf")
pdf(out_path, width = 6, height = 3.5)
for (sig_name in sig_names) {
  plot(cc_plts[[sig_name]])
}
dev.off()

#Heat map of medians

hm_data <- meta_stats %>% 
  rowwise() %>% 
  mutate(id = paste0(name, "__", Phase)) %>% 
  ungroup() %>% 
  select(id, cell_type_2, median) %>% 
  pivot_wider(names_from = id, values_from = median) %>% 
  column_to_rownames("cell_type_2") %>%
  t()

#invert columns order
hm_data <- hm_data[,ncol(hm_data):1]

rownames(hm_data) <- gsub("post_integration_pa_down_noRibo_", "PA-down_ribo-no", rownames(hm_data))
rownames(hm_data) <- gsub("post_integration_pa_down_", "PA-down_ribo-yes", rownames(hm_data))
rownames(hm_data) <- gsub("post_integration_pa_up_noRibo_", "PA-up_ribo-no", rownames(hm_data))
rownames(hm_data) <- gsub("post_integration_pa_up_", "PA-up_ribo-yes", rownames(hm_data))

rows_split <- gsub("_G1|_G2M|_S", "", rownames(hm_data))
cols_split <- groups %>% pull(cell_group)
cols_split <- factor(cols_split, levels = c("First_Line", "Lted", "TamR", "FulvR", "ER_Mutant"))

hm_data_sat <- hm_data
hm_data_sat[hm_data <= -0.25] <- -0.25
hm_data_sat[hm_data >= 0.25] <- 0.25

hm <- Heatmap(hm_data_sat,
              row_title_rot = 0,
              column_title_rot = 90,
              row_gap = unit(3, "mm"),
              heatmap_legend_param = list(title = "Module_Score"),
              row_names_max_width = unit(10, "cm"),
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              row_split = rows_split, 
              column_split = cols_split,
              border = TRUE,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8))

out_path <- paste0(out_folder, "/PA-scores.heatmap.pdf")
pdf(out_path, width = 7, height = 4)
draw(hm)
dev.off()

##Box plots using only cell_groups First_Line, TamR, Lted
##One box plots series per group

cgs <- c("First_Line", "TamR", "Lted")

cc_plts <- c()

for (cg in cgs) {
    
  meta_cg <- meta %>% dplyr::filter(cell_group == cg) %>% mutate(cell_type_2 = droplevels(cell_type_2))
  ids <- meta_cg %>% pull(cell_type_2) %>% levels()
  ids <- paste0("MCF7_", ids)

  cc_plts[["post_integration_pa_up"]][[cg]] <- ggplot(meta_cg, aes(x=cell_type_2, 
                                                                   y=post_integration_pa_up, 
                                                                   fill=cell_type_2)) +
    geom_boxplot() + 
    scale_fill_manual(values = cols_table[ids,1]) + 
    facet_wrap(~Phase) +
    coord_flip() +
    ylim(-0.75, 0.75) +
    ylab("post_integration_pa_up") +
    geom_hline(yintercept = 0, color = "black", lty=2) +
    theme_bw() + 
    theme(legend.position = "none")
  
  cc_plts[["post_integration_pa_down"]][[cg]] <- ggplot(meta_cg, aes(x=cell_type_2, 
                                                                   y=post_integration_pa_down, 
                                                                   fill=cell_type_2)) +
    geom_boxplot() + 
    scale_fill_manual(values = cols_table[ids,1]) + 
    facet_wrap(~Phase) +
    coord_flip() +
    ylim(-0.75, 0.75) +
    ylab("post_integration_pa_down") +
    geom_hline(yintercept = 0, color = "black", lty=2) +
    theme_bw() + 
    theme(legend.position = "none")
  
  cc_plts[["post_integration_pa_up_noRibo"]][[cg]] <- ggplot(meta_cg, aes(x=cell_type_2, 
                                                                          y=post_integration_pa_up_noRibo, 
                                                                          fill=cell_type_2)) +
    geom_boxplot() + 
    scale_fill_manual(values = cols_table[ids,1]) + 
    facet_wrap(~Phase) +
    coord_flip() +
    ylim(-0.75, 0.75) +
    ylab("post_integration_pa_up_noRibo") +
    geom_hline(yintercept = 0, color = "black", lty=2) +
    theme_bw() + 
    theme(legend.position = "none")
  
  cc_plts[["post_integration_pa_down_noRibo"]][[cg]] <- ggplot(meta_cg, aes(x=cell_type_2, 
                                                                            y=post_integration_pa_down_noRibo, 
                                                                            fill=cell_type_2)) +
    geom_boxplot() + 
    scale_fill_manual(values = cols_table[ids,1]) + 
    facet_wrap(~Phase) +
    coord_flip() +
    ylim(-0.75, 0.75) +
    ylab("post_integration_pa_down_noRibo") +
    geom_hline(yintercept = 0, color = "black", lty=2) +
    theme_bw() + 
    theme(legend.position = "none")
  
}

out_path <- paste0(out_folder, "/PA-scores.boxplots.grouped.ribo_yes.pdf")
pdf(out_path, width = 6, height = 9)
plot_grid(cc_plts$post_integration_pa_up$First_Line, 
          cc_plts$post_integration_pa_up$Lted,
          cc_plts$post_integration_pa_up$TamR, 
          cc_plts$post_integration_pa_down$First_Line, 
          cc_plts$post_integration_pa_down$Lted,
          cc_plts$post_integration_pa_down$TamR, 
          ncol = 1, 
          align="v")
dev.off()

out_path <- paste0(out_folder, "/PA-scores.boxplots.grouped.ribo_no.pdf")
pdf(out_path, width = 6, height = 9)
plot_grid(cc_plts$post_integration_pa_up_noRibo$First_Line, 
          cc_plts$post_integration_pa_up_noRibo$Lted,
          cc_plts$post_integration_pa_up_noRibo$TamR, 
          cc_plts$post_integration_pa_down_noRibo$First_Line, 
          cc_plts$post_integration_pa_down_noRibo$Lted,
          cc_plts$post_integration_pa_down_noRibo$TamR, 
          ncol = 1, 
          align="v")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~
#PA + Hallmarks Cor. Mat.
#########################

w1 <- grep("_pa_|HALLMARK", colnames(meta))
w2 <- grep("post_integration_", colnames(meta))
w <- intersect(w1, w2)
cm <- cor(meta[,w])
rownames(cm) <- gsub("HALLMARK_", "", rownames(cm))
colnames(cm) <- gsub("HALLMARK_", "", colnames(cm))
rownames(cm) <- gsub("post_integration_", "", colnames(cm))
colnames(cm) <- gsub("post_integration_", "", colnames(cm))

hm <- Heatmap(cm,
              row_title_rot = 0, 
              row_gap = unit(3, "mm"),
              heatmap_legend_param = list(title = "PCC"),
              row_names_max_width = unit(10, "cm"),
              border = TRUE,
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 8))

out_path <- paste0(out_folder, "/scores.cor_matrix.heatmap.pdf")
pdf(out_path, width = 11, height = 10)
draw(hm)
dev.off()

