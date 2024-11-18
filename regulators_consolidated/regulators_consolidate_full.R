## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: 
## regulators_consolidate_full.R
##
## Description: 
#~
## Consolidate regulators predictions
## from breast in vitro, breast in vivo
## prostate in vitro, prostate in vivo
## published DTP signatures
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
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)


#~~~~~~~~~~~~~~~~
#Prepare Datasets
#################

## Prostate in vitro (decoupleR; baseline subtracted)

d_pc <- read_tsv("../prostate_decipherseq_downstream/decoupleR.run_subtracted-paired/data.all.wide.anno.complete_info.txt")

#Exclude LNCaP95 + WM
#Basically one value first line, one value resistant, one value second line

d_pc_f <- d_pc %>% 
  dplyr::select(gene, `LNCaP_WM-7d`, LNCaP95, `LNCaP95_Enz-7d`)

## Breast in vitro (decoupleR; baseline subtracted)

d_bc <- read_tsv("../breast_regulators/decoupleR.stats.run_subtracted-paired/data.all.wide.anno.complete_info.txt")

#Use the aggregate first line, resistant, second line, mutant
#But also show the individual treatments (7d only for the second line)
#Already includes the analyses mentioned below on the LM neo-adjuvant

## Breast in vivo (pervasiveness; NES and p-value from fgsea)

d_bc_meta <- read_tsv("/Users/irosbarozzi/Documents/imperial/project_stress/neoadjuvant_202407/downstream.pervasiveness.results/fgsea.results.merge.txt")

#Dorothea
d_bc_meta <- d_bc_meta %>%
  dplyr::filter(gs_name == "Dorothea") %>%
  dplyr::select(-gs_name) %>%
  dplyr::rename(gene = pathway)


#~~~~~~~~~~~~~~~~~~~~
#Consolidate Datasets
#####################

#d_bc
#d_bc_meta
#d_pc_f

colnames(d_bc) <- paste0("BCa_", colnames(d_bc))
colnames(d_bc)[colnames(d_bc) == "BCa_gene"] <- "gene"

colnames(d_bc_meta) <- paste0("Perv_", colnames(d_bc_meta))
colnames(d_bc_meta)[colnames(d_bc_meta) == "Perv_gene"] <- "gene"

colnames(d_pc_f) <- paste0("PCa_", colnames(d_pc_f))
colnames(d_pc_f)[colnames(d_pc_f) == "PCa_gene"] <- "gene"

d <- d_bc %>% 
  left_join(d_bc_meta, by = "gene") %>% 
  left_join(d_pc_f, by = "gene")

## Statistics

# "BCa_First_line"  --> BCa_MCF7_first
# "BCa_Resistant"   --> BCa_MCF7_resistant
# "BCa_Second_line" --> BCa_MCF7_second
# "BCa_ESR1_mutant" --> BCa_MCF7_ERmutant

# "PCa_LNCaP_WM-7d"     --> PCa_LNCaP_first
# "PCa_LNCaP95"         --> PCa_LNCaP_resistant
# "PCa_LNCaP95_Enz-7d"  --> PCa_LNCaP_second

# Neo-adjuvant LM
# "BCa_neoadj_score" * -log10 "BCa_neoadj_pvalue_BH" --> BCa_NEO_new
# "BCa_neoadj_lasso_est"                             --> BCa_NEO_new_LASSO
# "BCa_neoadj_rf_imp"                                --> BCa_NEO_new_RF

# Pervasiveness
# "Perv_NEO_NES" * -log10 "Perv_NEO_padj" --> Perv_NEO
# "Perv_DTP_NES" * -log10 "Perv_DTP_padj" --> Perv_DTP

coi <-  c("Regulator", 
          "BCa_MCF7_first", "BCa_MCF7_resistant", "BCa_MCF7_second", "BCa_MCF7_ERmutant", 
          "PCa_LNCaP_first", "PCa_LNCaP_resistant", "PCa_LNCaP_second", 
          "BCa_NEO_new", "BCa_NEO_new_LASSO", "BCa_NEO_new_RF", 
          "Perv_NEO", "Perv_DTP")

d_stats <- d %>%
  rename(Regulator = gene) %>%
  rename(BCa_MCF7_first = BCa_First_line) %>%
  rename(BCa_MCF7_resistant = BCa_Resistant) %>%
  rename(BCa_MCF7_second = BCa_Second_line) %>%
  rename(BCa_MCF7_ERmutant = BCa_ESR1_mutant) %>%
  rename(PCa_LNCaP_first = `PCa_LNCaP_WM-7d`) %>%
  rename(PCa_LNCaP_resistant = PCa_LNCaP95) %>%
  rename(PCa_LNCaP_second = `PCa_LNCaP95_Enz-7d`) %>%
  mutate(BCa_NEO_new = BCa_neoadj_score * -log10(BCa_neoadj_pvalue_BH)) %>%
  rename(BCa_NEO_new_LASSO = BCa_neoadj_lasso_est) %>%
  rename(BCa_NEO_new_RF = BCa_neoadj_rf_imp) %>%
  mutate(Perv_NEO = Perv_NEO_NES * -log10(Perv_NEO_padj)) %>%
  mutate(Perv_DTP = Perv_DTP_NES * -log10(Perv_DTP_padj)) %>%
  dplyr::select(all_of(coi))

d_stats[is.na(d_stats)] <- 0

# Add flag for ML

d_stats <- d_stats %>% mutate(BCa_NEO_new_ML = BCa_NEO_new_LASSO != 0 | BCa_NEO_new_RF != 0)


#~~~~~~~~~~~~~
#Plots / Stats
##############

# Heat map of 
# BCa MCF7 + PCa LNCaP
# Annotate rows BCa_NEO + Pervasiveness

## Version 1 -- Top regulators in vitro (first line)

n_top <- 20

#using the same data, pre-filtered for significant ones
d_bc <- read_tsv("../breast_regulators/decoupleR_tests/decoupleR_tests.run_subtracted-paired/data.all.wide.anno.txt")
d_pc <- read_tsv("../prostate_decipherseq_downstream/fc_results/decoupleR.run_subtracted-paired/data.all.wide.anno.txt")

n_top_brca <- d_bc %>% arrange(-abs(First_line)) %>% head(n_top) %>% pull(gene)
n_top_prca <- d_pc %>% arrange(-abs(First_line)) %>% head(n_top) %>% pull(gene)

d_stats <- d_stats %>% mutate(Selection_1 = Regulator %in% c(n_top_brca, n_top_prca))

## Version 2 -- Top regulators in vivo (new NEO or Perv)

n_top_neo <- d_stats %>% arrange(-abs(BCa_NEO_new)) %>% head(n_top) %>% pull(Regulator)
n_top_perv <- d_stats %>% arrange(-abs(Perv_NEO)) %>% head(n_top) %>% pull(Regulator)

d_stats <- d_stats %>% mutate(Selection_2 = Regulator %in% c(n_top_neo, n_top_perv))

## Version 3 -- features selected in ML (neo)

d_stats <- d_stats %>% mutate(Selection_3 = BCa_NEO_new_ML)

## Heat maps

coi <-  c("BCa_MCF7_first", "BCa_MCF7_resistant", "BCa_MCF7_second", "BCa_MCF7_ERmutant", 
          "PCa_LNCaP_first", "PCa_LNCaP_resistant", "PCa_LNCaP_second")

selection_heights <- list()
selection_heights[["Selection_1"]] <- 9
selection_heights[["Selection_2"]] <- 11
selection_heights[["Selection_3"]] <- 19

for (col_name in c("Selection_1", "Selection_2", "Selection_3")) {
  
  d_stats_f <- d_stats %>% dplyr::filter(get(col_name))
  d_stats_mat <- d_stats_f %>% column_to_rownames(var = "Regulator") %>% select(all_of(coi))
  
  anno_ML <- d_stats_f$BCa_NEO_new_ML
  anno_NEO <- d_stats_f$BCa_NEO_new
  anno_NEO[anno_NEO >= 20] <- 20
  anno_NEO[anno_NEO <= -20] <- -20
  anno_Perv <- d_stats_f$Perv_NEO
  anno_Perv[anno_Perv >= 5] <- 5
  anno_Perv[anno_Perv <= -5] <- -5
  
  ha_cols <- list(BCa_NEO_new_ML = c("TRUE" = "orange", "FALSE" = "lightgrey"))
  
  ha <- rowAnnotation(BCa_NEO_new_ML = anno_ML,
                      BCa_NEO_new = anno_barplot(anno_NEO, ylim = c(-20,20)), 
                      Perv_NEO = anno_barplot(anno_Perv, ylim = c(-5,5)), 
                      annotation_name_gp= gpar(fontsize = 8),
                      gap = unit(1, "mm"),
                      col = ha_cols)
  
  csplit <- c(rep("Breast", 4), rep("Prostate", 3))
  
  set.seed(123)
  HM <- Heatmap(d_stats_mat,
                heatmap_legend_param = list(title = "TF-activity"),
                col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red")),
                column_title_rot = 90,
                row_title_rot = 0,
                clustering_distance_rows = "pearson",
                cluster_columns = FALSE,
                column_split = csplit, 
                border = TRUE,
                right_annotation = ha)
  
  out_file <- paste0("results.", col_name, ".heatmap.pdf")
  pdf(out_file, width = 5.5, height = selection_heights[[col_name]])
  draw(HM)
  dev.off()
  
}

## Version 4: either selection 1 (in vitro) or 2 (in vivo) + ESR1 & FOXA1 & AR
## Flagged by top in vitro bca, top in vitro pca, top in vivo, top pervasiveness

d_stats_f <- d_stats %>% dplyr::filter(Selection_1 | Selection_2 | Regulator %in% c("ESR1", "FOXA1", "AR"))
d_stats_f <- d_stats_f %>%
  mutate(Top_BrCa = Regulator %in% n_top_brca) %>%
  mutate(Top_PrCa = Regulator %in% n_top_prca) %>%
  mutate(Top_NeoNew = Regulator %in% n_top_neo) %>%
  mutate(Top_NeoPerv = Regulator %in% n_top_perv)

d_stats_mat <- d_stats_f %>% column_to_rownames(var = "Regulator") %>% select(all_of(coi))

anno_ML <- d_stats_f$BCa_NEO_new_ML
anno_NEO <- d_stats_f$BCa_NEO_new
anno_NEO[anno_NEO >= 20] <- 20
anno_NEO[anno_NEO <= -20] <- -20
anno_Perv <- d_stats_f$Perv_NEO
anno_Perv[anno_Perv >= 5] <- 5
anno_Perv[anno_Perv <= -5] <- -5

ha_cols <- list(Top_BrCa = c("TRUE" = "darkgreen", "FALSE" = "lightgrey"),
                Top_PrCa = c("TRUE" = "darkgreen", "FALSE" = "lightgrey"),
                Top_NeoNew = c("TRUE" = "purple", "FALSE" = "lightgrey"),
                Top_NeoPerv = c("TRUE" = "purple", "FALSE" = "lightgrey"),
                BCa_NEO_new_ML = c("TRUE" = "orange", "FALSE" = "lightgrey"))

ha <- rowAnnotation(Top_BrCa = d_stats_f$Top_BrCa,
                    Top_PrCa = d_stats_f$Top_PrCa,
                    Top_NeoNew = d_stats_f$Top_NeoNew,
                    BCa_NEO_new = anno_barplot(anno_NEO, ylim = c(-20,20)), 
                    Top_NeoPerv = d_stats_f$Top_NeoPerv,
                    Perv_NEO = anno_barplot(anno_Perv, ylim = c(-5,5)), 
                    BCa_NEO_new_ML = anno_ML,
                    annotation_name_gp= gpar(fontsize = 8),
                    gap = unit(1, "mm"),
                    col = ha_cols)

csplit <- c(rep("Breast", 4), rep("Prostate", 3))

set.seed(123)
HM <- Heatmap(d_stats_mat,
              heatmap_legend_param = list(title = "TF-activity"),
              col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red")),
              column_title_rot = 90,
              row_title_rot = 0,
              clustering_distance_rows = "pearson",
              cluster_columns = FALSE,
              column_split = csplit, 
              border = TRUE,
              right_annotation = ha)

out_file <- paste0("results.Selection_1+2.heatmap.pdf")
pdf(out_file, width = 6.5, height = 16)
draw(HM)
dev.off()

#same but with kmeans

set.seed(123)
HM <- Heatmap(d_stats_mat,
              heatmap_legend_param = list(title = "TF-activity"),
              col = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red")),
              column_title_rot = 90,
              row_title_rot = 0,
              row_km = 7, 
              row_km_repeats = 100,
              cluster_columns = FALSE,
              column_split = csplit, 
              border = TRUE,
              right_annotation = ha)

out_file <- paste0("results.Selection_1+2.heatmap.k=7.pdf")
pdf(out_file, width = 6.5, height = 15)
draw(HM)
dev.off()

## Potentially a Supplementary Table

write_tsv(x = d_stats, file = "results.txt")

## Borda count (~ sum of ranks)

d_stats_ranked <- d_stats %>% 
  dplyr::select(Regulator, BCa_MCF7_first, BCa_MCF7_second, PCa_LNCaP_first,
                BCa_NEO_new, Perv_NEO, Perv_DTP) %>%
  dplyr::mutate(BCa_MCF7_first_borda = rank(BCa_MCF7_first) - 1) %>%
  dplyr::mutate(BCa_MCF7_second_borda = rank(BCa_MCF7_second) - 1) %>%
  dplyr::mutate(PCa_LNCaP_first_borda = rank(PCa_LNCaP_first) - 1) %>%
  dplyr::mutate(BCa_NEO_new_borda = rank(BCa_NEO_new) - 1) %>%
  dplyr::mutate(Perv_NEO_borda = rank(Perv_NEO) - 1) %>%
  dplyr::mutate(Perv_DTP_borda = rank(Perv_DTP) - 1) %>% 
  dplyr::mutate(borda_count = BCa_MCF7_first_borda + BCa_MCF7_second_borda + PCa_LNCaP_first_borda + BCa_NEO_new_borda + Perv_NEO_borda + Perv_DTP_borda)

#maximum borda
borda_max <- (nrow(d_stats_ranked) - 1) * 6

#borda count over max
d_stats_ranked <- d_stats_ranked %>%
  dplyr::mutate(borda_fraction = borda_count / borda_max)

write_tsv(x = d_stats_ranked, file = "results.consolidated.txt")

#highlight top 10
reg_top10 <- d_stats_ranked %>% 
  dplyr::arrange(-borda_fraction) %>% 
  head(10) %>% 
  pull(Regulator)
d_stats_ranked <- d_stats_ranked %>%
  mutate(top = Regulator %in% reg_top10) %>%
  mutate(top_label = ifelse(top, Regulator, ""))

plt_scatter <- d_stats_ranked %>% 
  ggplot(aes(x = rank(borda_count), y = borda_fraction, col = top)) + 
  geom_point() + 
  theme_bw() + 
  xlab("Rank of Borda Count") + 
  ylab("Borda Count (over max)") +
  ylim(0,1) +
  theme(legend.position="none") +
  scale_colour_manual(values = c("darkgrey", "red"))

d_stats_ranked_sel <- d_stats_ranked %>% 
  dplyr::filter(top) %>% 
  mutate(Regulator = factor(Regulator, levels = rev(reg_top10)))

plt_bar <- d_stats_ranked_sel %>% 
  ggplot(aes(x = Regulator, xend = Regulator, y = 0.85, yend = borda_fraction)) + 
  geom_segment(col = "red", size = 5) +
  coord_flip() + 
  theme_bw() + 
  ylim(0.85,0.95) +
  ylab("Borda Count (over max)")

pdf("results.consolidated.top10.pdf", width = 5.5, height = 2.5)
grid.arrange(plt_scatter, plt_bar, widths = c(2, 1.8))
dev.off()

#second_line ~ first_line

cls <- rep("n.c.", nrow(d_stats_ranked))
cls[d_stats_ranked$BCa_MCF7_first >= 0.5 & d_stats_ranked$BCa_MCF7_second >= 0.5] <- "Up (both)"
cls[d_stats_ranked$BCa_MCF7_first <= -0.5 & d_stats_ranked$BCa_MCF7_second <= -0.5] <- "Down (both)"
cls[abs(d_stats_ranked$BCa_MCF7_first) < 0.5 & d_stats_ranked$BCa_MCF7_second >= 0.5] <- "Up (2nd only)"
cls[abs(d_stats_ranked$BCa_MCF7_first) < 0.5 & d_stats_ranked$BCa_MCF7_second <= -0.5] <- "Down (2nd only)"
d_stats_ranked <- d_stats_ranked %>% mutate(tf_class = cls)

scc <- cor.test(d_stats_ranked$BCa_MCF7_second, d_stats_ranked$BCa_MCF7_first, method = "spearman")

plt_title <- paste0("SCC = ", formatC(scc$estimate, 2), "; p = ", formatC(scc$p.value, 2))

plt <- d_stats_ranked %>% 
  mutate(gene_label = ifelse(cls != "n.c.", Regulator, NA)) %>% 
  ggplot(aes(x=BCa_MCF7_first, y=BCa_MCF7_second, color=tf_class, label=gene_label)) + 
  geom_hline(yintercept=0, linetype=2, color = "grey") +
  geom_vline(xintercept=0, linetype=2, color = "grey") +
  geom_point() +
  geom_text_repel(max.overlaps = 50) +
  scale_color_manual(values=c('darkblue', 'darkcyan', 'darkgrey', 'red', 'orange')) +
  xlim(c(-3.2, 3.2)) +
  ylim(c(-1.1, 1.1)) +
  ggtitle(plt_title) +
  theme_bw()

pdf("results.consolidated.brca-2nd_vs_brca-1st.pdf", width = 12, height = 8)
plot(plt)
dev.off()

