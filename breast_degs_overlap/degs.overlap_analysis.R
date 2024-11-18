## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: degs.overlap_analysis.R
##
## Description: 
#~
## Exploration of differential gene expression results for breast (and prostate) data. Emphasis on overlapping degs between different treatment strategies.
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

source("lib.R")

#~~~~~~~
# Params
########

secondLine_only_7d <- FALSE

#~~~~~~~~~~~~
#Prepare DEGs
#############

####~ 
####~ Load all markers from scRNA-seq experiments and add basic annotations
pa_up <- read_tsv("../data/PA_UP.txt", col_names = F) %>%
  dplyr::rename(gene = X1)
pa_down <- read_tsv("../data/PA_DOWN.txt", col_names = F) %>%
  dplyr::rename(gene = X1)

avg_expr_bca <- read_tsv("../breast_seurat/MCF7_cca_kit_5k_feat/avg_expr_breast.tsv")
avg_expr_pca <- read_tsv("../prostate_seurat/avg_expr_prostate.tsv")

markers_bca <- read_tsv("data/marker_genes.breast.txt")
deseq_bca <- readRDS("../breast_seurat/MCF7_cca_kit_5k_feat/deg_results_deseq.rds")
deseq_bca <- bind_rows(deseq_bca, .id = "id") %>%
  rename(symbol = feature, avg_logFC = log_fc, p_val_adj = padj, pct.1 = rate1, pct.2 = rate2) %>%
  select(-c(ave_expr, stat, pvalue, group1, group2)) %>%
  as_tibble() %>%
  filter(abs(avg_logFC) > 0.25) %>%
  filter(p_val_adj < 0.05) %>%
  mutate(PA_UP = ifelse(symbol %in% pa_up$gene, TRUE, FALSE),
         PA_DOWN = ifelse(symbol %in% pa_down$gene, TRUE, FALSE)) %>%
  dplyr::rename(avg_log2FC = avg_logFC) %>% 
  left_join(avg_expr_bca, by = c("symbol" = "gene"))


markers_pca <- read_tsv("data/marker_genes.prostate.txt")
deseq_pca <- readRDS("../prostate_seurat/deg_results_deseq.rds")
deseq_pca <- bind_rows(deseq_pca, .id = "id") %>%
  rename(symbol = feature, avg_logFC = log_fc, p_val_adj = padj, pct.1 = rate1, pct.2 = rate2) %>%
  select(-c(ave_expr, stat, pvalue, group1, group2)) %>%
  as_tibble() %>%
  filter(abs(avg_logFC) > 0.25) %>%
  filter(p_val_adj < 0.05) %>%
  mutate(PA_UP = ifelse(symbol %in% pa_up$gene, TRUE, FALSE),
         PA_DOWN = ifelse(symbol %in% pa_down$gene, TRUE, FALSE)) %>%
  dplyr::rename(avg_log2FC = avg_logFC) %>% 
  left_join(avg_expr_pca, by = c("symbol" = "gene"))

markers_bca <- markers_bca %>% dplyr::rename(avg_log2FC = avg_logFC) %>%
  mutate(avg_expr = NA)
markers_pca <- markers_pca %>% dplyr::rename(avg_log2FC = avg_logFC) %>%
  mutate(avg_expr = NA)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# inspect stats (old vs new)
############################
deseq_bca_n <- deseq_bca %>% mutate(type = "new", line = "breast")
deseq_pca_n <- deseq_pca %>% mutate(type = "new", line = "prostate")
markers_bca_o <- markers_bca %>% mutate(type = "old", line = "breast")
markers_pca_o <- markers_pca %>% mutate(type = "old", line = "prostate")

markers_o_n <- rbind(deseq_bca_n, deseq_pca_n, markers_bca_o, markers_pca_o) %>% 
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down"))

# number of DEGS
n_markers <- markers_o_n %>% group_by(type, line) %>% summarize(count = n())
n_markers_plt <- ggplot(n_markers, aes(x = type, y = count, fill = type)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  coord_flip() +
  facet_wrap(~line, scales = "free") +
  scale_fill_brewer(palette = "Accent") +
  theme(legend.position = "none")

n_markers_sample <- markers_o_n %>% group_by(type, line, id) %>% summarize(count = n())
n_markers_sample_plt <- ggplot(n_markers_sample, aes(x = id, y = count, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  facet_wrap(~line, scales = "free") +
  scale_fill_brewer(palette = "Accent") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0))

pdf("n_markers_old_new.pdf", width = 6, height = 2)
plot(n_markers_plt)
dev.off()

pdf("n_markers_sample_old_new.pdf", width = 12, height = 8)
plot(n_markers_sample_plt)
dev.off()

# FC and pval distribution
pval_distr <- ggplot(markers_o_n, aes(x = direction, y = -log10(p_val_adj), fill = type)) +
  geom_violin(position = position_dodge(0.8), trim = FALSE) +
  geom_boxplot(position = position_dodge(0.8), width = 0.05, outlier.shape = NA) +
  facet_wrap(~line) +
  scale_fill_brewer(palette = "Accent") +
  theme_classic()

fc_distr <- ggplot(markers_o_n, aes(x = direction, y = abs(avg_log2FC), fill = type)) +
  geom_violin(position = position_dodge(0.7), trim = FALSE) +
  geom_boxplot(position = position_dodge(0.7), width = 0.05, outlier.shape = NA) +
  facet_wrap(~line) +
  scale_fill_brewer(palette = "Accent") +
  theme_classic()

pdf("pval_markers_old_new.pdf", width = 7, height = 3)
plot(pval_distr)
dev.off()

pdf("fc_markers_old_new.pdf", width = 8, height = 4)
plot(fc_distr)
dev.off()

# average expression distribution
avg_expr_distr_log <- ggplot(markers_o_n %>% filter(type == "new"), aes(x = avg_expr, y = line, fill = line)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  scale_x_continuous(trans="log10") +
  scale_fill_brewer(palette = "Accent")

avg_expr_distr <- ggplot(markers_o_n %>% filter(type == "new"), aes(x = avg_expr, y = line, fill = line)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  scale_fill_brewer(palette = "Accent")

# binning average expression
expr_bins <- c("low" = 0.01, "medium" = 0.1, "high" = 1)
markers_o_n <- markers_o_n %>% 
  mutate(expr_bin = ifelse(avg_expr >= expr_bins["high"], "high", 
                           ifelse(avg_expr >= expr_bins["medium"], "medium", "low")))
overlap_tibbles <- list()
old_new_overlaps <- list()
for (cutoff in names(expr_bins)) {
  
  overlap_tibbles[[cutoff]] <- tibble(id = character(), direction = character(), type = character(), count = numeric(), avg_expr_threshold = numeric())
  for (id_use in unique(markers_o_n$id)) {
    m <- markers_o_n %>% filter(id == id_use)
    m_o <- m %>% filter(type == "old")
    m_n <- m %>% filter(type == "new")
    
    for (dir in c("up", "down")) {
      m_n_d <- m_n %>% filter(direction == dir, avg_expr >= expr_bins[cutoff])
      m_o_d <- m_o %>% filter(direction == dir)
      
      membership <- c("both" = length(intersect(m_n_d$symbol, m_o_d$symbol)),
                      "old_only" = length(setdiff(m_o_d$symbol, m_n_d$symbol)),
                      "new_only" = length(setdiff(m_n_d$symbol, m_o_d$symbol)))
      tib <- tibble("id" = id_use,
                    "direction" = dir,
                    "type" = names(membership), 
                    "count" = membership,
                    "avg_expr_threshold" = expr_bins[cutoff])
      
      overlap_tibbles[[cutoff]] <- rbind(overlap_tibbles[[cutoff]], tib)
    }
  }
  old_new_overlaps[[cutoff]] <- ggplot(overlap_tibbles[[cutoff]], aes(x = id, y = count, fill = type)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    scale_fill_brewer(palette = "Accent") +
    theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0)) +
    facet_wrap(~direction) +
    ggtitle(paste("Average expression cutoff:", expr_bins[cutoff]))
  
  write_tsv(x = overlap_tibbles[[cutoff]], file = paste0("deg_overlaps_old_new_", expr_bins[cutoff], ".tsv"))
}
write_tsv(x = markers_o_n, file = "markers_combined_old_new.tsv")

old_new_deg_overlaps_cp <- plot_grid(plotlist = old_new_overlaps, ncol = 1)
pdf("old_new_deg_overlaps.pdf", width = 10, height = 20)
plot(old_new_deg_overlaps_cp)
dev.off()


#~~~~~~~~~~~~~~~~~~~
# filter new markers
####################
#markers <- rbind(markers_bca, markers_pca)
#markers <- rbind(deseq_bca, deseq_pca)

markers <- markers_o_n %>% filter(type == "new", avg_expr >= expr_bins[["medium"]]) # only keep DEGs with an everage expression >= 0.1

ids_orig <- read_tsv("..//data/signatures_comparison_20210701.IDs.txt", col_names = FALSE) %>% pull(1)
ids_new <- paste("Novel_", sapply(ids_orig, function(x){strsplit(x, "_vs_")[[1]][1]}), sep = "")

for (i in 1:length(ids_orig)) {
  
  id_orig <- ids_orig[i]
  id_new <- ids_new[i]
  markers$id[markers$id == id_orig] <- id_new
  
}

markers <- markers %>% dplyr::filter(grepl("^Novel", id) & p_val_adj <= 0.05)

# Strip the Novel_* from the id(s)
markers$id <- gsub("Novel_", "", markers$id)

# Add column for regulation by therapy (up/down)
deg_col <- rep(NA, nrow(markers))
deg_col[markers$avg_log2FC < 0] <- "Down"
deg_col[markers$avg_log2FC > 0] <- "Up"
markers$therapy_regulation <- deg_col

# Add column PA (no/up/down)
pa_col <- rep("No", nrow(markers))
pa_col[markers$PA_DOWN] <- "Down"
pa_col[markers$PA_UP] <- "Up"
markers$PA_regulation <- pa_col

# Keep and arrange only the important info
markers <- markers %>% select(id, symbol, avg_log2FC, p_val_adj, pct.1, pct.2, therapy_regulation, PA_regulation)


####~ 
####~ Add further annotations

# Load samples metadata
pairs <- read_tsv("../data/pairs.txt")

# Add annotation PA / Adaptive

# First-line (easy, not PA is Adaptive)

gene_lists <- pairs %>% dplyr::filter(Line == "First") %>% pull(Sample) %>% unique()

markers_ext_first <- markers %>% dplyr::filter(id %in% gene_lists)
markers_ext_first <- markers_ext_first %>% mutate(gene_class = ifelse(PA_regulation == "No", "Adaptive", "PreAdaptive"))

#Second-line (get matched first-line, then stratify by Adaptive_FirstLine and Adaptive)

gene_lists <- pairs %>% dplyr::filter(Line == "Second") %>% pull(Sample) %>% unique()
markers_ext_second <- c()

for (gene_list in gene_lists) {
  
  matched_first <- pairs %>% dplyr::filter(Sample == gene_list) %>% pull(Matched_First)
  
  for (dir in c("Up", "Down")) {
    
    markers_first <- markers %>% dplyr::filter(id == matched_first & therapy_regulation == dir) %>% pull(symbol)
    markers_ext <- markers %>% dplyr::filter(id %in% gene_list & therapy_regulation == dir) %>% 
      rowwise() %>%
      mutate(gene_class = ifelse(PA_regulation == "No", ifelse(symbol %in% markers_first, "Adaptive_FirstLine", "Adaptive"), "PreAdaptive"))
    markers_ext_second  <- rbind(markers_ext_second, markers_ext)
    
  }
  
}

# Merge together
markers_ext <- rbind(markers_ext_first, markers_ext_second)

####~
####~ Filter conditions, if required
# Keep only 7d treatment for second line (and no mutant)

if (secondLine_only_7d) {
  
  pairs_first <- pairs %>% filter(Line == "First")
  pairs_second_filt <- pairs %>% 
    filter(Line == "Second") %>%
    filter(!grepl("2d", Sample))
  
  pairs <- rbind(pairs_first, pairs_second_filt)
  markers_ext <- markers_ext %>% filter(id %in% pairs$Sample)
  # and check if also need to filter markers (if used...)
  
}


####~ 
####~ Load other pre-calculated information about markers
# 
# d_all <- read_tsv("../pa_vs_acute.20220429_fig3_4SG/data/signatures.comparison.results.txt")
# 
# # Subset by columns of interests
# coi <- read_tsv("../pa_vs_acute.20220429_fig3_4SG/coi.txt", col_names = FALSE) %>% pull(1)
# d_all_coi <- d_all %>% dplyr::select(all_of(coi)) %>% dplyr::rename(symbol = Symbol)
# 
# # Strip the Novel_*
# colnames(d_all_coi) <- gsub("Novel_", "", colnames(d_all_coi))

####~ 
####~ Calculate fractions of markers across DEG comparisons (first-line and second-line)

markers_fractions <- list() 

for (line_sel in c("First", "Second")) {
  cols <- pairs %>% filter(Cell_Type == "Breast") %>% filter(Line == line_sel) %>% pull(Sample)
  col_name_up <- paste0("MCF7_UP_", line_sel, "Line_fraction")
  col_name_down <- paste0("MCF7_DOWN_", line_sel, "Line_fraction")
  
  markers_fractions[[line_sel]] <- 
    markers %>% 
    select(id, symbol, therapy_regulation) %>%                         # select relevant columns
    filter(id %in% cols) %>%                                           # filter by sample
    pivot_wider(names_from = id, values_from = therapy_regulation) %>% # pivot to wide format
    mutate(across(all_of(cols), ~replace_na(.x, "None"))) %>%          # replace NAs 
    rowwise() %>%
    mutate(total = length(cols),
           sum_up = sum(c_across(all_of(cols)) == "Up"),
           sum_down = sum(c_across(all_of(cols)) == "Down")) %>%       # calculate rowwise counts
    mutate(!!col_name_up := round(sum_up/total, 2),
           !!col_name_down := round(sum_down/total, 2)) %>%            # calculate fractions
    select(symbol, !!col_name_up, !!col_name_down)
}

# same as above, but only second line, stratified by treatment duration
line_sel <- "Second"
for (dur in c("2d", "7d")) {
  cols <- pairs %>% filter(Cell_Type == "Breast") %>% 
    mutate(duration = str_split_i(Sample, pattern = "-", -1)) %>% 
    filter(Line == line_sel, duration == dur) %>% pull(Sample)
  
  col_name_up <- paste0("MCF7_UP_", line_sel, "Line", dur, "_fraction")
  col_name_down <- paste0("MCF7_DOWN_", line_sel, "Line", dur, "_fraction")
  
  markers_fractions[[paste(line_sel, dur, sep = "_")]] <- 
    markers %>% 
    select(id, symbol, therapy_regulation) %>%                         # select relevant columns
    filter(id %in% cols) %>%                                           # filter by sample
    pivot_wider(names_from = id, values_from = therapy_regulation) %>% # pivot to wide format
    mutate(across(all_of(cols), ~replace_na(.x, "None"))) %>%          # replace NAs 
    rowwise() %>%
    mutate(total = length(cols),
           sum_up = sum(c_across(all_of(cols)) == "Up"),
           sum_down = sum(c_across(all_of(cols)) == "Down")) %>%       # calculate rowwise counts
    mutate(!!col_name_up := round(sum_up/total, 2),
           !!col_name_down := round(sum_down/total, 2)) %>%            # calculate fractions
    select(symbol, !!col_name_up, !!col_name_down)
}

####~ 
####~ ERa targets?

era_anno <- read_tsv("../data/er_target_classified.v2.txt")
era_anno_lst <- list()
era_anno_lst[["up"]] <- era_anno %>% dplyr::filter(ER_target_lenient == "up" & ER_target_ChIP) %>% pull(gene)
era_anno_lst[["down"]] <- era_anno %>% dplyr::filter(ER_target_lenient == "down" & ER_target_ChIP) %>% pull(gene)

for (anno_id in names(era_anno_lst)) {
  col_id <- paste("ERa_", anno_id, sep = "")
  markers_ext <- markers_ext %>% mutate(!!col_id := symbol %in% era_anno_lst[[anno_id]]) 
}
era_anno_col <- rep("not_a_target", rep(nrow(markers_ext)))
era_anno_col[markers_ext$ERa_up] <- "target_pos_reg"
era_anno_col[markers_ext$ERa_down] <- "target_neg_reg"
markers_ext <- markers_ext %>% mutate(ERa_annotation = era_anno_col)


####~ 
####~ Merge markers and calculated fraction information

markers_cons <- markers_ext %>%
  left_join(markers_fractions$First, by = "symbol") %>%
  left_join(markers_fractions$Second, by = "symbol") %>%
  left_join(markers_fractions$Second_2d, by = "symbol") %>%
  left_join(markers_fractions$Second_7d, by = "symbol")

write_tsv(markers_cons, file = "degs.annotated.txt")


#############

#~~~~~~~~~~~~~~~~~~~
#First-Line Analysis
####################

plots <- list()

gene_lists <- pairs %>% dplyr::filter(Line == "First") %>% pull(Sample) %>% unique()
# exclude LNCaP
gene_lists <- gene_lists[gene_lists != "LNCaP_WM-2d" & gene_lists != "LNCaP_WM-7d" & gene_lists != "MCF7-D538G_WM-2d"]

####~ 
####~ PreAdaptive vs Adaptive

plt_prefix <- "firstLine_PA"
pval_bp <- tibble()
for (dir in c("Up", "Down")) {
  
  stats <- markers_cons %>% 
    dplyr::filter(id %in% gene_lists & (therapy_regulation == dir)) %>% 
    select(id, gene_class) %>% 
    group_by(id, gene_class) %>%
    summarise(Genes = n()) %>%
    ungroup()
  
  for (plt_type in c("fill", "stack")) {
    
    plt <- ggplot(stats, aes(fill=gene_class, y=Genes, x=id)) + 
      geom_bar(position=plt_type, stat="identity", colour="black") +
      scale_fill_brewer(name = "", palette="Paired") +
      theme_bw() +
      coord_flip()
    
    plot_name <- paste(plt_prefix, "_", dir, "_barplot_", plt_type, sep = "")
    plots[[plot_name]] <- plt

  }
  
  #Fold-changes
  
  data_bp <- markers_cons %>% 
    dplyr::filter(id %in% gene_lists & (therapy_regulation == dir))
  
  plt <- ggplot(data_bp %>% as.data.frame(), mapping = aes(x = id, y = avg_log2FC, fill = gene_class)) + 
    geom_boxplot() +
    scale_fill_brewer(name = "", palette="Paired") +
    theme_bw() +
    coord_flip()
  
  for (test_id in (unique(data_bp$id))) {
    
    pa_fc <- data_bp %>% filter(id == test_id, gene_class == "PreAdaptive") %>% pull(avg_log2FC)
    a_fc <- data_bp %>% filter(id == test_id, gene_class == "Adaptive") %>% pull(avg_log2FC)
    test <- wilcox.test(x = pa_fc, y = a_fc, alternative = "two.sided", exact = FALSE)
    
    pval_bp <- rbind(
      pval_bp,
      tibble(test_id = paste("repr_vs_sel", test_id, dir, sep = "_"), pval = test$p.value, test = "Wilcoxon")
    )
    
  }
  
  plot_name <- paste(plt_prefix, "_", dir, "_fc_boxplot", sep = "")
  plots[[plot_name]] <- plt
  
}

####~ 
####~ PreAdaptive vs Adaptive - stratified by pervasiveness

plt_prefix <- "firstLine_PA_pervasiveness"

for (dir in c("Up", "Down")) {
  
  which_field <- ifelse(dir == "Up", "MCF7_UP_FirstLine_fraction", "MCF7_DOWN_FirstLine_fraction")
  
  markers_pre_filt <- markers_cons %>% 
    dplyr::filter(id %in% gene_lists & (therapy_regulation == dir))
  
  frac_levels <- markers_pre_filt %>% select(which_field) %>% pull() %>% unique() %>% sort()
  
  stats <- markers_pre_filt %>% 
    select(id, gene_class, which_field) %>% 
    dplyr::rename(Fraction = !!which_field) %>%
    mutate(Fraction = factor(Fraction, levels = frac_levels)) %>%
    group_by(id, gene_class, Fraction) %>%
    summarise(Genes = n()) %>%
    ungroup()
  
  for (plt_type in c("fill", "stack")) {
    
    plt <- ggplot(stats, aes(fill=Fraction, y=Genes, x=id)) + 
      geom_bar(position=plt_type, stat="identity", colour="black") +
      facet_wrap(~gene_class) +
      scale_fill_brewer(name = "", palette="Oranges") +
      theme_bw() +
      coord_flip()
    
    plot_name <- paste(plt_prefix, "_", dir, "_barplot_", plt_type, sep = "")
    plots[[plot_name]] <- plt

  }
  
  data_bp <- markers_cons %>% 
    dplyr::filter(id %in% gene_lists & (therapy_regulation == dir)) %>%
    dplyr::rename(Fraction = !!which_field) %>%
    mutate(Fraction = factor(Fraction, levels = frac_levels))
  
  plt <- ggplot(data_bp %>% as.data.frame(), mapping = aes(x = id, y = avg_log2FC, fill = Fraction)) + 
    geom_boxplot() +
    facet_wrap(~gene_class) +
    scale_fill_brewer(name = "", palette="Oranges") +
    theme_bw() +
    coord_flip()
  
  for (gene_class_sel in c("PreAdaptive", "Adaptive")) {
    
    for (id_sel in unique(data_bp$id)) {
      
      test_data <- data_bp %>% filter(id == id_sel, gene_class == gene_class_sel)
      test <- kruskal.test(x = test_data$avg_log2FC, test_data$Fraction)
      
      pval_bp <- rbind(
        pval_bp,
        tibble(test_id = paste("firstline_PA_pervasiveness", id_sel, gene_class_sel, dir, sep = "_"), pval = test$p.value, test = "Kruskal.wallis")
      )
      
    }
    
  }
  
  plot_name <- paste(plt_prefix, "_", dir, "_fc_boxplot", sep = "")
  plots[[plot_name]] <- plt

}

write_tsv(pval_bp, "overlaps.boxplots.tests.tsv")

####~ 
####~ Plots

plt_up_1 <- 
  cowplot::plot_grid(plots$firstLine_PA_Up_barplot_stack + ylim(0, 2500), 
                     plots$firstLine_PA_Up_barplot_fill, 
                     plots$firstLine_PA_Up_fc_boxplot,
                     labels = tolower(LETTERS)[1:3],
                     ncol = 3)

plt_up_2 <- 
  cowplot::plot_grid(plots$firstLine_PA_pervasiveness_Up_barplot_stack + ylim(0, 1700),
                     plots$firstLine_PA_pervasiveness_Up_barplot_fill,
                     plots$firstLine_PA_pervasiveness_Up_fc_boxplot,
                     labels = tolower(LETTERS)[4:6],
                     ncol = 3)

plt_down_1 <- 
  cowplot::plot_grid(plots$firstLine_PA_Down_barplot_stack + ylim(0, 2500), 
                     plots$firstLine_PA_Down_barplot_fill, 
                     plots$firstLine_PA_Down_fc_boxplot,
                     labels = tolower(LETTERS)[1:3],
                     ncol = 3)

plt_down_2 <- 
  cowplot::plot_grid(plots$firstLine_PA_pervasiveness_Down_barplot_stack + ylim(0, 1700),
                     plots$firstLine_PA_pervasiveness_Down_barplot_fill,
                     plots$firstLine_PA_pervasiveness_Down_fc_boxplot,
                     labels = tolower(LETTERS)[4:6],
                     ncol = 3)

plt_stats_f3a <- rbind(
  plots$firstLine_PA_Up_barplot_fill$data %>% mutate(dir = "up"),
  plots$firstLine_PA_Down_barplot_fill$data %>% mutate(dir = "down")
)

plt_stats_f3b <- rbind(
  plots$firstLine_PA_pervasiveness_Up_barplot_fill$data %>% mutate(dir = "up"),
  plots$firstLine_PA_pervasiveness_Down_barplot_fill$data %>% mutate(dir = "down")
)

write_tsv(plt_stats_f3a, "F3a_stats.tsv")
write_tsv(plt_stats_f3b, "F3b_stats.tsv")

pdf("overlaps.up.breast.first_line.stats.pdf", width = 13, height = 1.25)
plot(plt_up_1)
dev.off()

pdf("overlaps.up.breast.first_line.stats.pervasiveness.pdf", width = 16, height = 1.5)
plot(plt_up_2)
dev.off()

pdf("overlaps.down.breast.first_line.stats.pdf", width = 13, height = 1.25)
plot(plt_down_1)
dev.off()

pdf("overlaps.down.breast.first_line.stats.pervasiveness.pdf", width = 16, height = 1.5)
plot(plt_down_2)
dev.off()

####################

#~~~~~~~~~~~~~~~~~~~~
#Second-Line Analysis
#####################

gene_lists <- pairs %>% dplyr::filter(Line == "Second") %>% pull(Sample) %>% unique()
# exclude LNCaP95
gene_lists <- gene_lists[gene_lists != "LNCaP95_Enz-7d" & gene_lists != "LNCaP95_WM-7d"]

####~ 
####~ PreAdaptive vs Adaptive

plt_prefix <- "secondLine_PA"

for (dir in c("Up", "Down")) {
  
  stats <- markers_cons %>% 
    dplyr::filter(id %in% gene_lists & (therapy_regulation == dir)) %>% 
    select(id, gene_class) %>% 
    group_by(id, gene_class) %>%
    summarise(Genes = n()) %>%
    ungroup()
  
  for (plt_type in c("fill", "stack")) {
    
    plt <- ggplot(stats, aes(fill=gene_class, y=Genes, x=id)) + 
      geom_bar(position=plt_type, stat="identity", colour="black") +
      scale_fill_brewer(name = "", palette="Paired") +
      theme_bw() +
      coord_flip()
    
    plot_name <- paste(plt_prefix, "_", dir, "_barplot_", plt_type, sep = "")
    plots[[plot_name]] <- plt
    
  }
  
  #Fold-changes
  
  data_bp <- markers_cons %>% 
    dplyr::filter(id %in% gene_lists & (therapy_regulation == dir))
  
  plt <- ggplot(data_bp %>% as.data.frame(), mapping = aes(x = id, y = avg_log2FC, fill = gene_class)) + 
    geom_boxplot(outlier.stroke = 0.1) +
    scale_fill_brewer(name = "", palette="Paired") +
    theme_bw() +
    coord_flip()
  
  plot_name <- paste(plt_prefix, "_", dir, "_fc_boxplot", sep = "")
  plots[[plot_name]] <- plt
  
}

####~ 
####~ PreAdaptive vs Adaptive - stratified by pervasiveness

plt_prefix <- "secondLine_PA_pervasiveness"

for (dir in c("Up", "Down")) {
  
  which_field <- ifelse(dir == "Up", "MCF7_UP_SecondLine_fraction", "MCF7_DOWN_SecondLine_fraction")
  
  markers_pre_filt <-  markers_cons %>% 
    dplyr::filter(id %in% gene_lists & (therapy_regulation == dir))
  
  frac_levels <-markers_pre_filt  %>% select(which_field) %>% pull() %>% unique() %>% sort()
  
  stats <- markers_pre_filt %>% 
    select(id, gene_class, which_field) %>% 
    dplyr::rename(Fraction = !!which_field) %>%
    mutate(Fraction = factor(Fraction, levels = frac_levels)) %>%
    group_by(id, gene_class, Fraction) %>%
    summarise(Genes = n()) %>%
    ungroup()
  
  newPalette <- colorRampPalette(c("#fff5eb", "#fd8d3c", "#7f2704"))
  n_groups <- length(levels(stats$Fraction))
  
  for (plt_type in c("fill", "stack")) {
    
    plt <- ggplot(stats, aes(fill=Fraction, y=Genes, x=id)) + 
      geom_bar(position=plt_type, stat="identity", colour="black") +
      facet_wrap(~gene_class) +
      scale_fill_manual(values = newPalette(n_groups)) +
      theme_bw() +
      coord_flip()
    
    plot_name <- paste(plt_prefix, "_", dir, "_barplot_", plt_type, sep = "")
    plots[[plot_name]] <- plt
    
  }
  
  data_bp <- markers_cons %>% 
    dplyr::filter(id %in% gene_lists & (therapy_regulation == dir)) %>%
    dplyr::rename(Fraction = !!which_field) %>%
    mutate(Fraction = factor(Fraction, levels = frac_levels))
  
  plt <- ggplot(data_bp %>% as.data.frame(), mapping = aes(x = id, y = avg_log2FC, fill = Fraction)) + 
    geom_boxplot(outlier.stroke = 0.1) +
    facet_wrap(~gene_class) +
    scale_fill_manual(values = newPalette(length(unique(data_bp$Fraction)))) +
    theme_bw() +
    coord_flip()
  
  plot_name <- paste(plt_prefix, "_", dir, "_fc_boxplot", sep = "")
  plots[[plot_name]] <- plt
  
}

####~ 
####~ Plots

plt_up_1 <- 
  cowplot::plot_grid(plots$secondLine_PA_Up_barplot_stack, 
                     plots$secondLine_PA_Up_barplot_fill, 
                     plots$secondLine_PA_Up_fc_boxplot,
                     labels = tolower(LETTERS)[1:3],
                     ncol = 3)

plt_up_2 <- 
  cowplot::plot_grid(plots$secondLine_PA_pervasiveness_Up_barplot_stack,
                     plots$secondLine_PA_pervasiveness_Up_barplot_fill,
                     plots$secondLine_PA_pervasiveness_Up_fc_boxplot,
                     labels = tolower(LETTERS)[4:6],
                     ncol = 3)

plt_down_1 <- 
  cowplot::plot_grid(plots$secondLine_PA_Down_barplot_stack, 
                     plots$secondLine_PA_Down_barplot_fill, 
                     plots$secondLine_PA_Down_fc_boxplot,
                     labels = tolower(LETTERS)[1:3],
                     ncol = 3)

plt_down_2 <- 
  cowplot::plot_grid(plots$secondLine_PA_pervasiveness_Down_barplot_stack,
                     plots$secondLine_PA_pervasiveness_Down_barplot_fill,
                     plots$secondLine_PA_pervasiveness_Down_fc_boxplot,
                     labels = tolower(LETTERS)[4:6],
                     ncol = 3)

pdf("overlaps.up.breast.second_line.stats.pdf", width = 18, height = 4)
plot(plt_up_1)
dev.off()

pdf("overlaps.up.breast.second_line.stats.pervasiveness.pdf", width = 21, height = 4)
plot(plt_up_2)
dev.off()

pdf("overlaps.down.breast.second_line.stats.pdf", width = 18, height = 4)
plot(plt_down_1)
dev.off()

pdf("overlaps.down.breast.second_line.stats.pervasiveness.pdf", width = 21, height = 4)
plot(plt_down_2)
dev.off()


#####################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LTED vs. TAMR - pervasiveness
###############################

gene_lists <- pairs %>% dplyr::filter(Line == "Second") %>% pull(Sample) %>% unique()
# exclude LNCaP95
gene_lists <- gene_lists[gene_lists != "LNCaP95_Enz-7d" & gene_lists != "LNCaP95_WM-7d"]

plt_prefix <- "secondLine_treatment_pervasiveness"

# fraction levels for plotting
frac_levels <- c(markers_cons$MCF7_DOWN_SecondLine2d_fraction, markers_cons$MCF7_UP_SecondLine7d_fraction) %>% unique() %>% sort()

# # maximum number of genes to fix y-axis
# genes_all <- c()
# for (dir in c("Up", "Down")) {
#   which_field <- ifelse(dir == "Up", "MCF7_UP_SecondLine_fraction", "MCF7_DOWN_SecondLine_fraction")
#   genes <- markers_cons %>% 
#     select(id, !!which_field) %>%
#     dplyr::rename(Fraction = !!which_field) %>%
#     group_by(id, Fraction) %>%
#     summarize(count = n()) %>%
#     pull(count)
#   genes_all <- c(genes, genes_all)
# }
# genes_max <- max(genes_max)

plots_histograms <- list()
for (dir in c("Up", "Down")) {
  
  for (first_line in c("Lted", "TamR")) {
    
    second_lines <- if(first_line == "Lted") {
      second_lines <- c("Tam", "Fulv", "Palbo", "CDK7i")
    } else {
      second_lines <- c("WM", "Fulv", "Palbo", "CDK7i")
    }
    
    # different palettes for Lted and TamR?
    if (dir == "Up") {
      newPalette <- colorRampPalette(c("#fff5eb", "#fd8d3c", "#7f2704"))
      } else {
      newPalette <- colorRampPalette(c("#f0f7d4", "#66B032", "#092834"))
    }
    
    n_groups <- length(levels(stats$Fraction))
    
    for (second_line in second_lines) {
      
      if(dir == "Up") {
        which_fields <- c("MCF7_UP_SecondLine2d_fraction", "MCF7_UP_SecondLine7d_fraction")
      } else {
        which_fields <- c("MCF7_DOWN_SecondLine2d_fraction", "MCF7_DOWN_SecondLine7d_fraction")
      }
      
      gene_lists_filt <- gene_lists[grepl(pattern = first_line, x = gene_lists)]
      gene_lists_filt <- gene_lists_filt[grepl(pattern = second_line, x = gene_lists_filt)]
      
      markers_pre_filt <-  markers_cons %>% 
        dplyr::filter(id %in% gene_lists_filt & (therapy_regulation == dir))
      
      stats_2d <- markers_pre_filt %>%
        mutate(duration = str_split_i(id, pattern = "-", -1)) %>%
        filter(duration == "2d") %>%
        select(id, !!which_fields[1]) %>%
        dplyr::rename(Fraction = !!which_fields[1]) %>%
        mutate(Fraction = factor(Fraction, levels = frac_levels)) %>%
        group_by(id, Fraction) %>%
        summarise(Genes = n()) %>%
        ungroup()
      
      stats_7d <- markers_pre_filt %>%
        mutate(duration = str_split_i(id, pattern = "-", -1)) %>%
        filter(duration == "7d") %>%
        select(id, !!which_fields[2]) %>%
        dplyr::rename(Fraction = !!which_fields[2]) %>%
        mutate(Fraction = factor(Fraction, levels = frac_levels)) %>%
        group_by(id, Fraction) %>%
        summarise(Genes = n()) %>%
        ungroup()
      
      stats <- rbind(stats_2d, stats_7d) %>%
        mutate(duration = str_split_i(id, pattern = "-", -1))
        
      
      # stats <- markers_pre_filt %>% 
      #   select(id, !!which_fields[1], !!which_fields[2]) %>% 
      #   dplyr::rename(Fraction = !!which_field) %>%
      #   mutate(Fraction = factor(Fraction, levels = frac_levels)) %>%
      #   group_by(id, Fraction) %>%
      #   summarise(Genes = n()) %>%
      #   ungroup() %>%
      #   mutate(duration = str_split_i(id, pattern = "-", -1)) %>%
      #   mutate(duration = factor(duration, levels = c("2d", "7d")))
      
      plot_name <- paste(plt_prefix, dir, first_line, second_line, sep = "_")
      plots_histograms[[plot_name]] <- ggplot(stats, aes(x = Fraction, y = Genes, fill = duration)) +
        geom_bar(position = "identity", stat = "identity", colour = "black", alpha = 0.6) +
        #facet_wrap(~duration, dir = "v") +
        theme_bw() +
        #scale_fill_manual(values = newPalette(n_groups)) +
        #theme(legend.position = "none") +
        ggtitle(paste(first_line, second_line, sep = " ")) +
        xlab("Fraction") +
        ylab("Genes") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
      
    }
    
  }
  
}

plt_up <- cowplot::plot_grid(plots_histograms$secondLine_treatment_pervasiveness_Up_Lted_Tam,
                   plots_histograms$secondLine_treatment_pervasiveness_Up_TamR_WM,
                   plots_histograms$secondLine_treatment_pervasiveness_Up_Lted_Fulv,
                   plots_histograms$secondLine_treatment_pervasiveness_Up_TamR_Fulv,
                   plots_histograms$secondLine_treatment_pervasiveness_Up_Lted_Palbo,
                   plots_histograms$secondLine_treatment_pervasiveness_Up_TamR_Palbo,
                   plots_histograms$secondLine_treatment_pervasiveness_Up_Lted_CDK7i,
                   plots_histograms$secondLine_treatment_pervasiveness_Up_TamR_CDK7i,
                   labels = tolower(LETTERS),
                   ncol = 2)

plt_down <- cowplot::plot_grid(plots_histograms$secondLine_treatment_pervasiveness_Down_Lted_Tam,
                   plots_histograms$secondLine_treatment_pervasiveness_Down_TamR_WM,
                   plots_histograms$secondLine_treatment_pervasiveness_Down_Lted_Fulv,
                   plots_histograms$secondLine_treatment_pervasiveness_Down_TamR_Fulv,
                   plots_histograms$secondLine_treatment_pervasiveness_Down_Lted_Palbo,
                   plots_histograms$secondLine_treatment_pervasiveness_Down_TamR_Palbo,
                   plots_histograms$secondLine_treatment_pervasiveness_Down_Lted_CDK7i,
                   plots_histograms$secondLine_treatment_pervasiveness_Down_TamR_CDK7i,
                   labels = tolower(LETTERS),
                   ncol = 2)


pdf("overlaps.up.breast.second_line.histograms.pervasiveness.pdf", width = 5, height = 7)
plot(plt_up)
dev.off()

pdf("overlaps.down.breast.second_line.histograms.pervasiveness.pdf", width = 5, height = 7)
plot(plt_down)
dev.off()

#~~~
# summary lineplot (7d only)
colors <- read_tsv("../data/colors_for_umap.txt")

line_plot <- list()
plot_data <- list()
for (dir in c("Up", "Down")) {
  which_field <- ifelse(dir == "Up", "MCF7_UP_SecondLine7d_fraction", "MCF7_DOWN_SecondLine7d_fraction")
  
  stats <- markers_cons %>% 
    dplyr::filter(id %in% gene_lists & (therapy_regulation == dir)) %>%
    select(id, symbol, !!which_field) %>%
    dplyr::rename(Fraction = !!which_field) %>%
    mutate(duration = str_split_i(id, pattern = "-", -1)) %>%
    filter(duration == "7d") %>%
    group_by(id, Fraction) %>%
    summarize(Genes = n()) %>%
    group_by(id) %>%
    mutate(id = str_split(id, pattern = "_", n = 2)[[1]][2]) %>%
    mutate(id = str_split_i(id, pattern = "-", i = 1)) %>%
    mutate(id = factor(id, levels = c("Lted_Tam", "Lted_Fulv", "Lted_Palbo", "Lted_CDK7i",
                                      "TamR_WM", "TamR_Fulv", "TamR_Palbo", "TamR_CDK7i"))) %>%
    mutate(label = if_else(Fraction == 0.12, as.character(id), NA_character_)) %>% 
    mutate(fraction_Genes = Genes/sum(Genes))
  
  col_code <- colors %>%
    group_by(description) %>%
    mutate(description = str_split(description, pattern = "_", n = 2)[[1]][2]) %>%
    mutate(description = str_split_i(description, pattern = "-", i = 1)) %>%
    filter(description %in% unique(stats$id))
  
  col_vec <- col_code$color
  names(col_vec) <- col_code$description
  
  line_plot[[dir]] <- ggplot(stats, aes(x = Fraction, y = Genes, color = id)) +
    geom_line(linewidth = 1.2) +
    theme_bw() +
    scale_color_manual(values = col_vec) +
    geom_label_repel(aes(label = label),
                     size = 3,
                     fontface = "bold",
                     nudge_x = -0.2,
                     segment.linetype = "dashed",
                     segment.alpha = .7,
                     segment.curvature = -0.1,
                     segment.ncp = 3,
                     segment.angle = 20,
                     direction = "y",
                     box.padding = 0.12,
                     hjust = 0.25,
                     na.rm = TRUE) +
    theme(legend.position = "none") +
    ggtitle(paste("Overlapping fractions", dir, "(7d treatment)"))
  
  line_plot[[paste0(dir, "_normalized")]] <- ggplot(stats, aes(x = Fraction, y = fraction_Genes, color = id)) +
    geom_line(linewidth = 1.2) +
    theme_bw() +
    scale_color_manual(values = col_vec) +
    geom_label_repel(aes(label = label),
                     size = 3,
                     fontface = "bold",
                     nudge_x = -0.2,
                     segment.linetype = "dashed",
                     segment.alpha = .7,
                     segment.curvature = -0.1,
                     segment.ncp = 3,
                     segment.angle = 20,
                     direction = "y",
                     box.padding = 0.12,
                     hjust = 0.25,
                     na.rm = TRUE) +
    theme(legend.position = "none") +
    ggtitle(paste("Overlapping fractions", dir, "(7d treatment)"))
  plot_data[[dir]] <- stats %>% mutate(dir = dir)
}

line_grid <- cowplot::plot_grid(line_plot$Up,
                                line_plot$Down,
                                line_plot$Up_normalized,
                                line_plot$Down_normalized,
                                ncol = 2,
                                labels = tolower(LETTERS))

plt_data_f3fg <- rbind(plot_data$Up, plot_data$Down) %>% select(-label)
write_tsv(plt_data_f3fg, "plt_data_f3fg.tsv")

pdf("overlaps.breast.second_line_7d.lineplots.pdf", width = 12, height = 8)
plot(line_grid)
dev.off()

#####################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Discordant Regulation (PA vs Adaptive)
#######################################

stats_disc <- markers_cons %>% 
  dplyr::filter(PA_regulation != "No") %>% 
  select(id, therapy_regulation, PA_regulation) %>% 
  mutate(cls = paste(PA_regulation, therapy_regulation, sep = "_")) %>% 
  select(id, cls) %>% 
  table() %>% 
  as_tibble()

plt_stack <- ggplot(stats_disc, aes(fill=cls, y=n, x=id)) + 
  geom_bar(position="stack", stat="identity", colour="black") +
  scale_fill_brewer(name = "", palette="Paired") +
  theme_bw() +
  coord_flip()

plt_fill <- ggplot(stats_disc, aes(fill=cls, y=n, x=id)) + 
  geom_bar(position="fill", stat="identity", colour="black") +
  scale_fill_brewer(name = "", palette="Paired") +
  theme_bw() +
  coord_flip()

plt_merge <- 
  cowplot::plot_grid(plt_stack,
                     plt_fill,
                     labels = tolower(LETTERS)[1:2],
                     ncol = 2)

pdf("overlaps.discordant.PA_vs_Adaptive.pdf", width = 11, height = 5)
plot(plt_merge)
dev.off()


#######################################

#~~~~~~~~~~~~~~~~~~~~~~~~
#Prostate Cancer Analysis
#########################

#Show it later in a different way


#########################

#~~~~~~~~~~~~~~~~
#PCA of the genes
#################

markers_pca <- markers_cons

#-log10P
markers_pca <- markers_pca %>% mutate(p_val_adj = -log10(p_val_adj)) #%>%
  
# added this to replace infinite values with maximum value of that column
maximum <- max(markers_pca$p_val_adj[is.finite(markers_pca$p_val_adj)])
markers_pca <- markers_pca %>%
  mutate(p_val_adj = ifelse(is.infinite(p_val_adj), maximum, p_val_adj)) 

#avg_log2FC pivot
fc_pivot <- markers_pca %>% select(symbol, id, avg_log2FC) %>% pivot_wider(names_from = id, values_from = avg_log2FC)
fc_pivot[is.na(fc_pivot)] <- 0

#avg_log2FC pivot
padj_pivot <- markers_pca %>% select(symbol, id, p_val_adj) %>% pivot_wider(names_from = id, values_from = p_val_adj)
padj_pivot[is.na(padj_pivot)] <- 1

#join back to pervasiveness
markers_pca_perv <- markers_pca %>% select(symbol, grep("_fraction", colnames(markers_pca))) %>% unique() %>%
  replace(is.na(.), 0) # added this to replace NAs
markers_pca <- markers_pca_perv %>% full_join(fc_pivot, by = "symbol")
markers_pca <- markers_pca %>% left_join(padj_pivot, by = "symbol", suffix = c("_log2fc", "_padj"))

####~ 
####~ Prepare data for PCA

normalizeMinMax <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

#Based on fold-changes, p-values, pervasiveness

#scaling all using normalizeMinMax
#log2FC & -log10P: scale 0-1
#pervasiveness: should stay as it is

for (cname in colnames(markers_pca)) {
  if (cname != "symbol") {
    markers_pca <- markers_pca %>% mutate(!!cname := normalizeMinMax(get(cname)))
  }
}


####~ 
####~ PCA

res <- prcomp(markers_pca %>% column_to_rownames("symbol"), center = FALSE)

#% of variance explained

res_sdev <- tibble(PC = colnames(res$rotation), sdev_frac = res$sdev / sum(res$sdev))
write_tsv(res_sdev, file = "overlaps.pca.sdev.txt")

#also consolidate the results with genes metadata [rotation + genes metadata]
markers_anno <- markers_pca %>% left_join(markers_cons %>% select(symbol, therapy_regulation, PA_regulation, gene_class, ERa_annotation) %>% unique(), by = "symbol")

# markers_anno <- markers_cons %>% 
#   select(symbol, therapy_regulation, PA_regulation, gene_class, ERa_annotation,
#          MCF7_UP_FirstLine_fraction, MCF7_UP_SecondLine_fraction, MCF7_DOWN_FirstLine_fraction, MCF7_DOWN_SecondLine_fraction) %>% 
#   arrange(symbol) %>% 
#   unique()

res_anno <- as_tibble(res$x)
res_anno <- res_anno %>% mutate(symbol = markers_pca %>% pull(symbol))
res_anno <- res_anno %>% left_join(markers_anno, by = "symbol")
res_anno <- res_anno %>% select(symbol, therapy_regulation, PA_regulation, gene_class, ERa_annotation,
                                MCF7_UP_FirstLine_fraction, MCF7_UP_SecondLine_fraction, 
                                MCF7_DOWN_FirstLine_fraction, MCF7_DOWN_SecondLine_fraction, everything())
write_tsv(res_anno, file = "overlaps.pca.results.txt")

p_a <- ggplot(res_anno, aes(x=PC1, y=PC2)) + 
  geom_point(aes(col=MCF7_UP_FirstLine_fraction)) + 
  facet_wrap(~PA_regulation)
p_b <- ggplot(res_anno, aes(x=PC1, y=PC2)) + 
  geom_point(aes(col=MCF7_UP_SecondLine_fraction)) + 
  facet_wrap(~PA_regulation)
p_c <- ggplot(res_anno, aes(x=PC1, y=PC2)) + 
  geom_point(aes(col=MCF7_DOWN_FirstLine_fraction)) + 
  facet_wrap(~PA_regulation)
p_d <- ggplot(res_anno, aes(x=PC1, y=PC2)) + 
  geom_point(aes(col=MCF7_DOWN_SecondLine_fraction)) + 
  facet_wrap(~PA_regulation)

plt_merge <- cowplot::plot_grid(
  p_a, p_b, p_c, p_d,
  labels = tolower(LETTERS)[1:4],
  ncol = 1)

pdf("overlaps.pca.results.pc2_vs_pc1.pdf", width = 8, height = 12)
plot(plt_merge)
dev.off()

plt_data_f3c <- res_anno %>% select(symbol, 
                                    MCF7_UP_FirstLine_fraction, 
                                    MCF7_UP_SecondLine_fraction, 
                                    MCF7_DOWN_FirstLine_fraction, 
                                    MCF7_DOWN_SecondLine_fraction) %>%
  distinct()

write_tsv(plt_data_f3c, "plt_data_f3c.tsv")

p_a <- ggplot(plt_data_f3c, aes(x=MCF7_UP_FirstLine_fraction %>% as.factor, y=MCF7_UP_SecondLine_fraction)) + 
  geom_boxplot() +
  theme_bw() +
  ylim(0, 1) +
  xlab("First-line pervasiveness") +
  ylab("Second-line pervasiveness") +
  ggtitle("Up-regulated")

p_b <- ggplot(plt_data_f3c, aes(x=MCF7_DOWN_FirstLine_fraction %>% as.factor, y=MCF7_DOWN_SecondLine_fraction)) + 
  geom_boxplot() +
  theme_bw() +
  ylim(0, 1) +
  xlab("First-line pervasiveness") +
  ylab("Second-line pervasiveness") +
  ggtitle("Down-regulated")

plt_merge <- cowplot::plot_grid(
  p_a, p_b,
  labels = tolower(LETTERS)[1:4],
  ncol = 1)

pdf("overlaps.pca.results.pervasiveness_boxplot.pdf", width = 3, height = 6)
plot(plt_merge)
dev.off()

#################
