
gs <- list()


#~~~~~~~~~~~~~~~~
#MSigDB Hallmarks
#################

#gsea sets
gsea_H <- msigdbr(species = "Homo sapiens", category = "H")
gsea_H_t2g <- gsea_H %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
gsea_TF <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD")
gsea_TF_t2g <- gsea_TF %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

for (gs_name in unique(gsea_H_t2g$gs_name)) {
  w <- gsea_H_t2g$gs_name == gs_name
  gs[[gs_name]] <- gsea_H_t2g[w,"gene_symbol"]
}


#~~~~~~~~~~~
#Pre-Adapted
############

in_f <- "../data/PA.DOWN.txt"
gs[["PA_DOWN"]] <- read_tsv(in_f, col_names = FALSE) %>% pull(1) %>% sort()

in_f <- "../data/PA.UP.txt"
gs[["PA_UP"]] <- read_tsv(in_f, col_names = FALSE) %>% pull(1) %>% sort()

#without cell-cycle (CC) genes

cc.genes <- readRDS("../data/cc.genes_seurat.rds")

cc_genes_merge <- unique(c(unlist(cc.genes), gs$HALLMARK_G2M_CHECKPOINT))

w <- ! gs[["PA_DOWN"]] %in% cc_genes_merge
gs[["PA_DOWN_NO_CC"]] <- gs[["PA_DOWN"]][w]

#Almost the same, only 2 genes removed
#w <- ! gs[["PA_UP"]] %in% cc_genes_merge
#gs[["PA_UP_NO_CC"]] <- gs[["PA_UP"]][w]

#also formatted for GSEA
gsea_pa_t2g <- c()
for (name in c("PA_UP", "PA_DOWN", "PA_DOWN_NO_CC")) {
  gsea_pa_t2g <- rbind(gsea_pa_t2g, cbind(name, gs[[name]]))
}
gsea_pa_t2g <- data.frame(gs_name = gsea_pa_t2g[,1],
                          gene_symbol = gsea_pa_t2g[,2],
                          stringsAsFactors = FALSE)

#PA (tibble)
gs_tib <- list()
gs_tib[["pa"]] <- gsea_pa_t2g %>% as_tibble()

#PA (list)
gs_lst <- list()
for (set_id in names(gs_tib)) {
  gs_lst[[set_id]] <- list()
  ids <- gs_tib[[set_id]] %>% pull(gs_name) %>% unique()
  for (id in ids) {
    gs_lst[[set_id]][[id]] <- gs_tib[[set_id]] %>% filter(gs_name == id) %>% pull(gene_symbol)
  }
}
