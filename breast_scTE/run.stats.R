## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: run.stats.R
##
## Description: 
#~
## Investigate transposable elements
## (TE) expression in different
## conditions
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

library(Seurat)
library(SeuratData)
library(SeuratDisk)

anno_rmsk <- read_tsv("rmsk.ids.txt")
anno_rmsk_ids <- anno_rmsk %>% pull(Subfamily)


#~~~~~~~~~
#Load Data
##########

prefixes <- system("ls *.h5seurat", intern = TRUE)
prefixes <- gsub(".h5seurat", "", prefixes)

so_ind <- list()

for (prefix in prefixes) {
  inFile <- paste0(prefix, ".h5seurat")
  ID <- gsub("scTE_out_", "", prefix)
  so_ind[[ID]] <- LoadH5Seurat(inFile, meta.data = TRUE)
  cell_ids <- colnames(so_ind[[ID]]$RNA)
  so_ind[[ID]]@meta.data <- data.frame(orig.ident = rep(ID, length(cell_ids)))
  rownames(so_ind[[ID]]@meta.data) <- cell_ids
}

so_full <- merge(so_ind[["MCF7_RM"]], y = so_ind[["MCF7_WM_2d"]], project = "TEST")
for (ID in names(so_ind)) {
  if (! ID %in% c("MCF7_RM", "MCF7_WM_2d")) {
    so_full <- merge(so_full, y = so_ind[[ID]], project = "TEST")
  }
}

# Sub sampling
so <- subset(x = so_full, downsample = 10000)

# Re-order the identities
so@meta.data$orig.ident <- factor(so@meta.data$orig.ident, levels = c("MCF7_RM", "MCF7_WM_2d", "LTED", "MCF7_TAM_2d", "MCF7-TAMR", "MCF7-D538G_RM_2d", "MCF7-D538G_WM_2d"))


#~~~~~~~~~~~~~~~
#Quality Control
################

## Add Metadata

w_rmsk <- rownames(so@assays$RNA) %in% anno_rmsk_ids

#Counts total
so <- AddMetaData(
  object = so,
  metadata = colSums(so@assays$RNA),
  col.name = 'nCount_RNA'
)
#Counts genes
so <- AddMetaData(
  object = so,
  metadata = colSums(so@assays$RNA[!w_rmsk,]),
  col.name = 'nCount_RNA_genes'
)
#Counts repeats
so <- AddMetaData(
  object = so,
  metadata = colSums(so@assays$RNA[w_rmsk,]),
  col.name = 'nCount_RNA_repeats'
)

#Features total
so <- AddMetaData(
  object = so,
  metadata = colSums(so@assays$RNA[,] > 0),
  col.name = 'nFeature_RNA'
)
#Features genes
so <- AddMetaData(
  object = so,
  metadata = colSums(so@assays$RNA[!w_rmsk,] > 0),
  col.name = 'nFeature_RNA_genes'
)
#Features repeats
so <- AddMetaData(
  object = so,
  metadata = colSums(so@assays$RNA[w_rmsk,] > 0),
  col.name = 'nFeature_RNA_repeats'
)

#Mito percentage
so <- AddMetaData(
  object = so,
  metadata = PercentageFeatureSet(so, pattern = "^MT-"),
  col.name = 'percent.mt'
)

#VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", 
#                         "nCount_RNA_genes", "nCount_RNA_repeats", 
#                         "nFeature_RNA_genes", "nFeature_RNA_repeats"),
#        ncol = 7,
#        pt.size = 0,
#        group.by = "orig.ident")

#Filter

nFeature_RNA_min <- 250   # min number of genes expressed per cell
nCount_RNA_min <- 500     # min number of transcripts expressed per cell
nCount_RNA_max <- 20000   # max number of transcripts expressed per cell
percent.mt_max <- 15      # max % of mito reads

so <- subset(so, percent.mt <= percent.mt_max &
                 nFeature_RNA >= nFeature_RNA_min &
                 nCount_RNA >= nCount_RNA_min &
                 nCount_RNA <= nCount_RNA_max)

## Add fraction of fragments from repetitive elements

so@meta.data$nFeature_RNA_repeats_frac <- so@meta.data$nFeature_RNA_repeats / so@meta.data$nFeature_RNA

p <- VlnPlot(so, 
             features = c("nFeature_RNA_repeats_frac"), 
             group.by = "orig.ident", 
             pt.size = 0, 
             y.max = 0.25)

pdf("results.pdf", width = 6, height = 4)
plot(p)
dev.off()

## Exclude mito genes and set min.features to 100

keep <- grep("^MT-", rownames(so@assays$RNA), invert = TRUE)

so.filt <- CreateSeuratObject(counts = so@assays$RNA[keep,], project = "TEST", meta.data = so@meta.data, min.features = nFeature_RNA_min)


#~~~~~~~~~~~~~~~
#Seurat Pipeline
################

so.filt <- NormalizeData(so.filt, normalization.method = "LogNormalize", scale.factor = 5000)
so.filt <- FindVariableFeatures(so.filt, selection.method = "vst", nfeatures = 3000)
so.filt <- ScaleData(so.filt)
so.filt <- RunPCA(so.filt, npcs = 100)

Idents(so.filt) <- so.filt@meta.data$orig.ident
ids <- Idents(so.filt) %>% levels()
ids <- ids[ids != "MCF7_RM"]
all.markers.list <- list()
for (id_i in ids) {
  all.markers.list[[id_i]] <- FindMarkers(so.filt, ident.1 = id_i, ident.2 = "MCF7_RM")
}

#Collapse
all.markers <- c()
for (id_i in ids) {
  id_tib <- all.markers.list[[id_i]] %>% 
    rownames_to_column(var = "gene") %>% 
    as_tibble() %>% 
    mutate(cluster = id_i)
  all.markers <- rbind(all.markers, id_tib)
}

#Add flag & annotation for repeat
all.markers <- all.markers %>% 
  mutate(repeat_flag = gene %in% anno_rmsk_ids) %>%
  left_join(anno_rmsk %>% dplyr::rename(gene = Subfamily), by = "gene")

write_tsv(all.markers, file = "results.markers.txt")

#Annotate rmsk.ids.txt with all.markers

clusters <- c("MCF7_WM_2d", "LTED", "MCF7_TAM_2d", "MCF7-TAMR", "MCF7-D538G_RM_2d", "MCF7-D538G_WM_2d")

anno_rmsk_up <- anno_rmsk

for (i  in 1:length(clusters)) {
  
  cluster_i <- clusters[i]
  
  cluster_i_mrk <- all.markers %>% 
    dplyr::filter(cluster == cluster_i & repeat_flag & p_val_adj <= 0.05 & avg_log2FC > 0) %>% 
    pull(gene)
  
  id <- paste0("up_", cluster_i)
  
  anno_rmsk_up <- anno_rmsk_up %>% 
    mutate(!!id := Subfamily %in% cluster_i_mrk)

}

#save full table
write_tsv(anno_rmsk_up, file = "results.rmsk.ids.markers.txt")

#save the contingency table for L1
anno_rmsk_up_L1_cont <- anno_rmsk_up %>% 
  dplyr::filter(Family == "L1") %>% 
  select(up_MCF7_WM_2d, `up_MCF7-D538G_RM_2d`) %>% 
  table() %>% 
  as_tibble()
write_tsv(anno_rmsk_up_L1_cont, file = "results.rmsk.ids.markers.L1.CT.txt")


#~~~~~~~~~~~~~~~~~~~~~~~~~
#Repeat Families - Markers
##########################

repeats_stats <- all.markers %>% 
  mutate(sign = ifelse(avg_log2FC > 0, "up", "down")) %>%
  dplyr::filter(cluster != "MCF7_RM" & p_val_adj <= 0.05 & repeat_flag) %>% 
  group_by(Family, cluster, sign) %>% 
  summarise(n_sig = n()) %>%
  ungroup()

repeats_stats$cluster <- factor(repeats_stats$cluster, 
                                levels = c("MCF7_WM_2d", "LTED", "MCF7_TAM_2d", "MCF7-TAMR", "MCF7-D538G_RM_2d", "MCF7-D538G_WM_2d"))

feats <- c("Alu", "L1", "ERV1", "ERVK", "ERVL")

for (lvl in levels(repeats_stats$cluster)) {
  for (feat in feats) {
    for (sign_i in c("up", "down")) {
      w_row <- repeats_stats %>% dplyr::filter(cluster == lvl & Family == feat & sign == sign_i)
      if (nrow(w_row) == 0) {
        repeats_stats <- rbind(repeats_stats, c(feat, lvl, sign_i, 0))
      }
    }
  }
}

repeats_stats$n_sig <- as.numeric(repeats_stats$n_sig)

p <- ggplot(data=repeats_stats %>% dplyr::filter(Family %in% feats), aes(x=Family, y=n_sig, fill=cluster)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values = c("gold", "red", "pink", "magenta", "lightgreen", "darkgreen")) +
  facet_wrap(~sign) +
  coord_flip() +
  theme_bw()

pdf("results.markers.rm_family.barplots.pdf", width = 6.5, height = 3)
plot(p)
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~
#Repeat Families - Collated
###########################

#Create a copy of the Seurat object
so_cp <- so.filt

for (slot_name in c("data", "scale.data")) {
  
  #Extract the counts for the classes and families of repeats
  if (slot_name == "data") {
    em <- so_cp@assays$RNA@data
  } else if ((slot_name == "scale.data")) {
    em <- so_cp@assays$RNA@scale.data
  }
  w <- rownames(em) %in% anno_rmsk$Subfamily
  em_pvt <- em[w,] %>% 
    as.data.frame() %>% 
    rownames_to_column("gene") %>% 
    as_tibble() %>% 
    pivot_longer(cols = c(-"gene"), names_to = "cell")
  em_pvt <- em_pvt %>% left_join(anno_rmsk, by = c("gene" = "Subfamily"))
  
  #Summarize
  em_pvt_fam <- em_pvt %>% 
    group_by(Family, cell) %>% 
    summarise(expr = sum(value))
  
  #Add the new summarized expression values into the Seurat Object
  em_fam <- em_pvt_fam %>% 
    ungroup() %>% 
    pivot_wider(names_from = cell, values_from = expr)
  em_fam_df <- em_fam %>% 
    select(-c("Family")) %>% 
    as.data.frame()
  rownames(em_fam_df) <- em_fam %>% 
    pull(Family)
  
  if (slot_name == "data") {
    so_cp@assays$RNA@data <- rbind(so_cp@assays$RNA@data, as.matrix(em_fam_df))
  } else if (slot_name == "scale.data") {
    so_cp@assays$RNA@scale.data <- rbind(so_cp@assays$RNA@scale.data, as.matrix(em_fam_df))
  }
  
}

#Re-run markers (only on those, using "features")
all.markers.rm_fam <- FindAllMarkers(so_cp, features = em_fam %>% pull(Family), min.pct	= 0.01, logfc.threshold = 0.1)
write_tsv(all.markers.rm_fam, file = "results.markers.rm_family.txt")

#Plot results for each family

feats <- c("Alu", "L1", "ERV1", "ERVK", "ERVL")
p_vln_feats_id <- list()

for (feat in feats) {
  
  p_vln_feats_id[[feat]] <- VlnPlot(
    so_cp, features = feat, pt.size = 0, group.by = "orig.ident") +
    stat_summary(fun = median, geom='point', size = 10, colour = "black", shape = 95) +
    theme(legend.position = 'none') +
    theme(plot.title = element_text(size=12)) + 
    theme(plot.margin = unit(c(1,1.5,1,1.5), "cm")) + 
    theme(axis.text.x = element_text(size=10)
    )
  
}

pdf("results.rm_family.violin.id.pdf", width = 3, height = 4)
for (feat in feats) {
  plot(p_vln_feats_id[[feat]])
}
dev.off()

