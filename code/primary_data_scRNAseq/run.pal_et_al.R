
library("tidyverse")
library("Seurat")
library("msigdbr")
library("RColorBrewer")

vioplot_ylim <- c(-0.3,0.9)
vioplot_genes_ylim <- c(-1,6)

#Colors for plotting clusters
cols <- list()
cols[["clusters"]] <- list()
#Basal, LP, ML, Unknown
cols[["clusters"]][["NormEpiSub"]] <- c("darkgreen", "blue", "orange", "darkred")
#Basal, LP, ML, Unknown
cols[["clusters"]][["ERTotalTum"]] <- c("blue", "orange", "darkgreen")
#CAFs, Immune, T1, T2
cols[["clusters"]][["Male"]] <- c("darkred", "darkgreen", "blue", "orange")

#~~~~~~~~~~~~~~~~
#Data Preparation
#################

in_folder <- "./"

#SeuratObject_NormEpiSub.rds
#SeuratObject_ERTotalTum.rds
#SeuratObject_Male.rds
#can be downloaded from figshare
#https://figshare.com/articles/dataset/Data_R_code_and_output_Seurat_Objects_for_single_cell_RNA-seq_analysis_of_human_breast_tissues/17058077

dois <- c("NormEpiSub", "ERTotalTum", "Male")

data <- list()

for (doi in dois) {
  in_F <- paste(in_folder, "SeuratObject_", doi, ".rds", sep = "")
  data[[doi]] <- readRDS(in_F)
}

## Load Metadata

in_F <- paste(in_folder, "data_pal_et_al/embj2020107333-sup-0006-tableev4.txt", sep = "")
metadata <- read_tsv(in_F)
metadata <- metadata %>% 
  mutate(group = `Sample Name`) %>% 
  mutate(group = gsub("-", "_", group))

## Prepare signatures of interest

#Hallmark gene sets

gs_hallmark <- list()

gsea_H <- msigdbr(species = "Homo sapiens", category = "H")
gsea_H_t2g <- gsea_H %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
for (gs_name in unique(gsea_H_t2g$gs_name)) {
  w <- gsea_H_t2g$gs_name == gs_name
  gs_hallmark[[gs_name]] <- gsea_H_t2g[w,"gene_symbol"]
}

#PA

gs_pa <- list()

in_f <- "../data/PA_UP.txt"
gs_pa[["PA_UP"]] <- read_tsv(in_f, col_names = FALSE) %>% pull(1) %>% sort()

in_f <- "../data/PA_DOWN.txt"
gs_pa[["PA_DOWN"]] <- read_tsv(in_f, col_names = FALSE) %>% pull(1) %>% sort()

#PA down without cell-cycle (CC) genes
cc_genes_merge <- unique(c(unlist(cc.genes), gs_hallmark$HALLMARK_G2M_CHECKPOINT))
w <- ! gs_pa[["PA_DOWN"]] %in% cc_genes_merge
gs_pa[["PA_DOWN_NO_CC"]] <- gs_pa[["PA_DOWN"]][w]

## Add scores

#PA signature

for (doi in dois) {
  for (sig in names(gs_pa)) {
    data[[doi]] <- AddModuleScore(object = data[[doi]],
                                  features = list(gs_pa[[sig]]),
                                  ctrl = 10,
                                  name = sig,
                                  seed = 1234)
    w <- colnames(data[[doi]]@meta.data) == paste(sig, "1", sep ="")
    colnames(data[[doi]]@meta.data)[w] <- sig
  }
}

## Rename clusters according to the paper annotation

#NormEpiSub

n <- nrow(data[["NormEpiSub"]]@meta.data)
data[["NormEpiSub"]]@meta.data$cluster <- rep("LP", n)
w <- data[["NormEpiSub"]]@meta.data$seurat_clusters == 1
data[["NormEpiSub"]]@meta.data$cluster[w] <- "ML"
w <- data[["NormEpiSub"]]@meta.data$seurat_clusters == 2
data[["NormEpiSub"]]@meta.data$cluster[w] <- "Basal"
w <- data[["NormEpiSub"]]@meta.data$seurat_clusters == 3
data[["NormEpiSub"]]@meta.data$cluster[w] <- "Unknown"

#ERTotalTum

n <- nrow(data[["ERTotalTum"]]@meta.data)
data[["ERTotalTum"]]@meta.data$cluster <- rep("1", n)
w <- data[["ERTotalTum"]]@meta.data$seurat_clusters == 1
data[["ERTotalTum"]]@meta.data$cluster[w] <- "2"
w <- data[["ERTotalTum"]]@meta.data$seurat_clusters == 2
data[["ERTotalTum"]]@meta.data$cluster[w] <- "3"

#Male

n <- nrow(data[["Male"]]@meta.data)
data[["Male"]]@meta.data$cluster <- rep("Tumour_1", n)
w <- data[["Male"]]@meta.data$seurat_clusters == 1
data[["Male"]]@meta.data$cluster[w] <- "Tumour_2"
w <- data[["Male"]]@meta.data$seurat_clusters == 2
data[["Male"]]@meta.data$cluster[w] <- "Immune"
w <- data[["Male"]]@meta.data$seurat_clusters == 3
data[["Male"]]@meta.data$cluster[w] <- "CAFs"

## Add metadata

#NormEpiSub; add pre-/post- menopause

meta_add <- as_tibble(data[["NormEpiSub"]]@meta.data) %>% 
  select(group) %>% 
  left_join(metadata, by = "group") %>%
  select(Menopause, Parity)

data[["NormEpiSub"]]@meta.data$Menopause <- meta_add$Menopause
data[["NormEpiSub"]]@meta.data$Parity <- meta_add$Parity

#################

#~~~~~~~~~~~~~~~~~~~~~~~
#Normal Epithelial Cells
########################

#Basal, LP, ML, Unknown
clusters_cols <- cols[["clusters"]][["NormEpiSub"]]
#Post, Pre
menopause_cols <- c("darkred", "white")
#Nulli, Parous
parity_cols <- c("white", "darkred")

#UMAP
res <- DimPlot(data[["NormEpiSub"]], group.by = "cluster", cols = clusters_cols)
pdf("NormEpiSub.UMAP.pdf", width=7, height=6)
plot(res)
dev.off()

#Violin plots by cluster
pdf("NormEpiSub.PA.vioplots.cluster.pdf", width=4, height=6)
for (sig in names(gs_pa)) {
   res <- VlnPlot(data[["NormEpiSub"]], features = c(sig), pt.size = 0, group.by = "cluster", cols = clusters_cols) + 
     stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
     theme(legend.position = 'bottom') +
     theme(plot.title = element_text(size=10)) + 
     theme(plot.margin = unit(c(1,1.5,1,1.5), "cm")) + 
     theme(axis.text.x = element_text(size=10)) +
     scale_y_continuous(limits = vioplot_ylim)
   plot(res)
}
dev.off()

#Violin plots by cluster + menopause
pdf("NormEpiSub.PA.vioplots.cluster_and_menopause.pdf", width=6, height=6)
for (sig in names(gs_pa)) {
  res <- VlnPlot(data[["NormEpiSub"]], features = c(sig), pt.size = 0, group.by = "cluster", split.by = "Menopause", cols = menopause_cols) + 
    theme(legend.position = 'bottom') +
    theme(plot.title = element_text(size=10)) + 
    theme(plot.margin = unit(c(1,1.5,1,1.5), "cm")) + 
    theme(axis.text.x = element_text(size=10)) +
    scale_y_continuous(limits = vioplot_ylim)
  plot(res)
}
dev.off()

#Violin plots by cluster + parity
pdf("NormEpiSub.PA.vioplots.cluster_and_parity.pdf", width=6, height=6)
for (sig in names(gs_pa)) {
  res <- VlnPlot(data[["NormEpiSub"]], features = c(sig), pt.size = 0, group.by = "cluster", split.by = "Parity", cols = parity_cols) + 
    theme(legend.position = 'bottom') +
    theme(plot.title = element_text(size=10)) + 
    theme(plot.margin = unit(c(1,1.5,1,1.5), "cm")) + 
    theme(axis.text.x = element_text(size=10)) +
    scale_y_continuous(limits = vioplot_ylim)
  plot(res)
}
dev.off()

#Metadata for scatterplots
p_data <- data[["NormEpiSub"]]@meta.data %>% as_tibble()

#Scatter plots PA UP ~ DOWN
pdf("NormEpiSub.PA.scatterplot.pdf", width=6, height=4.5)
p <- ggplot(p_data, aes(y=PA_UP, x=PA_DOWN, color=cluster)) +
  geom_point(size = 0.3) + 
  theme_bw() +
  theme(plot.title = element_text(size=10)) +
  xlim(vioplot_ylim) +
  ylim(vioplot_ylim) +
  scale_color_manual(values=clusters_cols)
plot(p)
dev.off()

#Scatter plots PA UP ~ DOWN: split by clusters
pdf("NormEpiSub.PA.scatterplot.clusters.pdf", width=7, height=5.5)
p <- ggplot(p_data, aes(y=PA_UP, x=PA_DOWN, color=cluster)) +
  facet_wrap("cluster") +
  geom_point(size = 0.5) + 
  theme_bw() +
  theme(plot.title = element_text(size=10)) +
  xlim(vioplot_ylim) +
  ylim(vioplot_ylim) +
  stat_smooth(color = "black", method = "lm", se = TRUE, fullrange = TRUE) +
  stat_cor(label.x.npc = "center", 
           label.y.npc = "top",
           size = 3) +
  scale_color_manual(values=clusters_cols)
plot(p)
dev.off()

#FeaturePlot(data[["ERTotalTum"]], features = c("PA_UP"), min.cutoff = -1, max.cutoff = 2)

########################

#~~~~~~~~~~~~
#Tumour Cells
#############

#Number of cells per donor per cluster
ERTotalTum_stats <- data[["ERTotalTum"]]@meta.data %>% 
  select(group, cluster) %>% 
  table() %>% 
  as_tibble()
write_tsv(ERTotalTum_stats, file = "ERTotalTum.stats.txt")

#Basal, LP, ML, Unknown
clusters_cols <- cols[["clusters"]][["ERTotalTum"]]

#UMAP
res <- DimPlot(data[["ERTotalTum"]], group.by = "cluster", cols = clusters_cols)
pdf("ERTotalTum.UMAP.pdf", width=7, height=6)
plot(res)
dev.off()

#Violin plots by cluster
pdf("ERTotalTum.PA.vioplots.cluster.pdf", width=3, height=6)
for (sig in names(gs_pa)) {
  res <- VlnPlot(data[["ERTotalTum"]], features = c(sig), pt.size = 0, group.by = "cluster", cols = clusters_cols) + 
    stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
    theme(legend.position = 'bottom') +
    theme(plot.title = element_text(size=10)) + 
    theme(plot.margin = unit(c(1,1.5,1,1.5), "cm")) + 
    theme(axis.text.x = element_text(size=10)) +
    scale_y_continuous(limits = vioplot_ylim)
  plot(res)
}
dev.off()

#Metadata for scatterplots
p_data <- data[["ERTotalTum"]]@meta.data %>% as_tibble()

#Scatter plots PA UP ~ DOWN
pdf("ERTotalTum.PA.scatterplot.pdf", width=6, height=4.5)
p <- ggplot(p_data, aes(y=PA_UP, x=PA_DOWN, color=cluster)) +
  geom_point(size = 0.3) + 
  theme_bw() +
  theme(plot.title = element_text(size=10)) +
  scale_color_manual(values=clusters_cols)
plot(p)
dev.off()

#Scatter plots PA UP ~ DOWN: split by clusters
pdf("ERTotalTum.PA.scatterplot.clusters.pdf", width=9.5, height=3)
p <- ggplot(p_data, aes(y=PA_UP, x=PA_DOWN, color=cluster)) +
  facet_wrap("cluster") +
  geom_point(size = 0.5) + 
  theme_bw() +
  theme(plot.title = element_text(size=10)) +
  stat_smooth(color = "black", method = "lm", se = TRUE, fullrange = TRUE) +
  stat_cor(label.x.npc = "center", 
           label.y.npc = "top",
           size = 3) +
  scale_color_manual(values=clusters_cols)
plot(p)
dev.off()

#FeaturePlot(data[["NormEpiSub"]], features = c("PA_UP"), min.cutoff = -1, max.cutoff = 2)

#############

#~~~~~~~~~~~~~~~~~~~
#+ Male Tumour Cells
####################

#CAFs, Immune, T1, T2
clusters_cols <- cols[["clusters"]][["Male"]]

#UMAP
res <- DimPlot(data[["Male"]], group.by = "cluster", cols = clusters_cols)
pdf("Male.UMAP.pdf", width=7, height=6)
plot(res)
dev.off()

#Violin plots by cluster
pdf("Male.PA.vioplots.cluster.pdf", width=4, height=6)
for (sig in names(gs_pa)) {
  res <- VlnPlot(data[["Male"]], features = c(sig), pt.size = 0, group.by = "cluster", cols = clusters_cols) + 
    stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
    theme(legend.position = 'bottom') +
    theme(plot.title = element_text(size=10)) + 
    theme(plot.margin = unit(c(1,1.5,1,1.5), "cm")) + 
    theme(axis.text.x = element_text(size=10)) +
    scale_y_continuous(limits = vioplot_ylim)
  plot(res)
}
dev.off()

#Metadata for scatterplots
p_data <- data[["Male"]]@meta.data %>% as_tibble() %>% filter(cluster %in% c("Tumour_1", "Tumour_2"))

#Scatter plots PA UP ~ DOWN
pdf("Male.PA.scatterplot.pdf", width=6, height=4.5)
p <- ggplot(p_data, aes(y=PA_UP, x=PA_DOWN, color=cluster)) +
  geom_point(size = 0.3) + 
  theme_bw() +
  theme(plot.title = element_text(size=10)) +
  xlim(vioplot_ylim) +
  ylim(vioplot_ylim) +
  scale_color_manual(values=clusters_cols[3:4])
plot(p)
dev.off()

#Scatter plots PA UP ~ DOWN: split by clusters
pdf("Male.PA.scatterplot.clusters.pdf", width=7, height=5.5)
p <- ggplot(p_data, aes(y=PA_UP, x=PA_DOWN, color=cluster)) +
  facet_wrap("cluster") +
  geom_point(size = 0.5) + 
  theme_bw() +
  theme(plot.title = element_text(size=10)) +
  xlim(vioplot_ylim) +
  ylim(vioplot_ylim) +
  stat_smooth(color = "black", method = "lm", se = TRUE, fullrange = TRUE) +
  stat_cor(label.x.npc = "center", 
           label.y.npc = "top",
           size = 3) +
  scale_color_manual(values=clusters_cols[3:4])
plot(p)
dev.off()

#FeaturePlot(data[["Male"]], features = c("PA_UP"), min.cutoff = -1, max.cutoff = 2)

####################

#~~~~~~~~~~~~~~~~~~
#Housekeeping Genes
###################

hk_genes <- c("GAPDH", "RPL26", "RPL36")

for (doi in dois) {
  out_F <- paste(doi, ".HK.vioplots.pdf", sep = "")
  pdf(out_F, width=4, height=6)
  for (i in 1:length(hk_genes)) {
    res <- VlnPlot(data[[doi]], features = hk_genes[i], pt.size = 0, group.by = "cluster", cols = cols[["clusters"]][[doi]]) + 
      stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
      theme(legend.position = 'bottom') +
      theme(plot.title = element_text(size=10)) + 
      theme(plot.margin = unit(c(1,1.5,1,1.5), "cm")) + 
      theme(axis.text.x = element_text(size=10)) +
      scale_y_continuous(limits = vioplot_genes_ylim)
    plot(res)
  }
  dev.off()
}

###################

#~~~~~~~~~~~~~~~~~~~~~~~~
#Marker genes of interest
#########################

feats <- c("MKI67", "ESR1", "FOXA1", "SOX4", "KLF6", "CD44", "TFF1")

for (doi in dois) {
  out_F <- paste(doi, ".goi.vioplots.pdf", sep = "")
  pdf(out_F, width=4, height=6)
  for (i in 1:length(feats)) {
    res <- VlnPlot(data[[doi]], features = feats[i], pt.size = 0, group.by = "cluster", cols = cols[["clusters"]][[doi]]) + 
      stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
      theme(legend.position = 'bottom') +
      theme(plot.title = element_text(size=10)) + 
      theme(plot.margin = unit(c(1,1.5,1,1.5), "cm")) + 
      theme(axis.text.x = element_text(size=10))
    plot(res)
  }
  dev.off()
}

#########################

