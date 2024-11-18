
source("downstream.libs.R")

outFolder <- "downstream.neo.treatment.results/"
cmd <- paste("mkdir", outFolder)
system(cmd)

#~~~~~~~~~~~~~~~~
#Data Preparation
#################

## Hallmark Gene Sets

gs_hallmark <- list()

gsea_H <- msigdbr(species = "Homo sapiens", category = "H")
gsea_H_t2g <- gsea_H %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
for (gs_name in unique(gsea_H_t2g$gs_name)) {
  w <- gsea_H_t2g$gs_name == gs_name
  gs_hallmark[[gs_name]] <- gsea_H_t2g[w,"gene_symbol"]
}

## PA Gene Sets

gs_pa <- list()

in_f <- "../data/PA_UP.txt"
gs_pa[["PA_UP"]] <- read_tsv(in_f, col_names = FALSE) %>% pull(1) %>% sort()

in_f <- "../data/PA_DOWN.txt"
gs_pa[["PA_DOWN"]] <- read_tsv(in_f, col_names = FALSE) %>% pull(1) %>% sort()

#PA down without cell-cycle (CC) genes
cc.genes <- readRDS("../data/cc.genes_seurat.rds")
cc_genes_merge <- unique(c(unlist(cc.genes), gs_hallmark$HALLMARK_G2M_CHECKPOINT))
w <- ! gs_pa[["PA_DOWN"]] %in% cc_genes_merge
gs_pa[["PA_DOWN_NO_CC"]] <- gs_pa[["PA_DOWN"]][w]

#ours
inF <- "sleuth.results/sleuth.results.txt"
sleuth_table <- read_tsv(file = inF)
neo_ranked_genes <- sleuth_table %>% 
  mutate(score = -log10(qval)*sign(b)) %>%
  arrange(score)
neo_ranked_genes_lst <- neo_ranked_genes %>% pull(score)
names(neo_ranked_genes_lst) <- neo_ranked_genes %>% pull(target_id)

#published
source("downstream.published.dataPrep.R")


#~~~~
#GSEA
#####

#Our Neo-adjuvant vs other sets

gs_all_lst <- list()
gs_all_lst[["pa"]] <- gs_pa
gs_all_lst[["hallmark"]] <- gs_hallmark
gs_all_lst[["neo_published"]] <- gs_neo

for (set_id in names(gs_all_lst)) {
  
  fgseaRes <- list()
  fgseaRes[[set_id]] <- fgsea(pathways = gs_all_lst[[set_id]],
                              stats    = neo_ranked_genes_lst,
                              eps      = 0.0,
                              minSize  = 10,
                              maxSize  = 1500)
  
  outF <- paste(outFolder, "/neo.gsea.", set_id, ".txt", sep="")
  fwrite(fgseaRes[[set_id]], file = outF, sep="\t", sep2=c("", ",", ""))
  
  #Plotting the significant ones
  
  sets <- fgseaRes[[set_id]] %>% filter(padj <= 0.05) %>% pull(pathway)
  
  outF <- paste(outFolder,"/neo.gsea.", set_id, ".significant.pdf", sep = "")
  pdf(outF, width = 5, height = 3)
  
  for (set in sets) {
    NES <- fgseaRes[[set_id]] %>% filter(pathway == set) %>% pull(NES) %>% round(2)
    plot_title <- paste(set, "\nNES = ", NES, sep = "")
    res <- plotEnrichment(gs_all_lst[[set_id]][[set]], neo_ranked_genes_lst) + 
      labs(title=plot_title) +
      theme(plot.title = element_text(size=8))
    plot(res)
  }
  
  dev.off()
  
}


#~~~~~~~~~~~~~~~~~~~~~~~
#Neo-adjuvant comparison
########################

#Bin genes

input_lst <- -neo_ranked_genes_lst  #invert so it can be overlaid on GSEA

br <- quantile(input_lst, seq(0,1,0.05))
bin <- as.numeric(cut(input_lst, breaks = br, include.lowest = TRUE))

neo_binned <- tibble(gene = names(input_lst),
                     val = input_lst,
                     bin = bin)

for (set_id in names(gs_neo)) {
  neo_binned <- neo_binned %>% 
    mutate(!!set_id := gene %in% gs_neo[[set_id]])
}

#Compute overlap per dataset per bin
neo_binned_anno <- neo_binned %>%
  group_by(bin) %>%
  summarise_at(names(gs_neo), function(x,tot){x = sum(x)})

#Normalize the overlap by the corresponding total number of genes in each signature
rnames <- neo_binned_anno %>% pull(bin)
neo_binned_anno_df <- neo_binned_anno %>% select(-c("bin")) %>% as.data.frame()
rownames(neo_binned_anno_df) <- rnames
neo_binned_anno_df <- t(neo_binned_anno_df) / sapply(gs_neo, length)

#Heat maps

#Sort the data first
var_sort_1 <- c("Neoadj")
var_sort_2 <- c("UP", "DOWN")
w <- c()
row_groups <- c()
for (v1 in var_sort_1) {
  for (v2 in var_sort_2) {
    w1 <- grep(v1, rownames(neo_binned_anno_df))
    w2 <- grep(v2, rownames(neo_binned_anno_df))
    wi <- intersect(w1, w2)
    w <- c(w, wi)
    lbl <- paste(v1, "_", v2, sep = "")
    row_groups <- c(row_groups, rep(lbl, length(wi)))
  }
}

hm_data <- neo_binned_anno_df[w,]

HM <- Heatmap(hm_data,
              heatmap_legend_param = list(title = "Fraction of Genes"),
              col = colorRamp2(c(0, 0.25, 0.5), c("white", "orange", "red")),
              row_split = row_groups,
              row_gap = unit(c(3), "mm"),
              row_names_gp = grid::gpar(fontsize = 10),
              column_names_gp = grid::gpar(fontsize = 0),
              cluster_columns = FALSE,
              cluster_rows = FALSE,
              cluster_row_slices = FALSE,
              row_title_rot = 0, 
              column_title_rot = 90,
              border = TRUE)

#move the caption
outF <- paste(outFolder,"/neo.overlap.neo_published.heatmap.pdf", sep = "")
pdf(outF, width=7, height=1.5)
draw(HM)
dev.off()

