# Adapted by Stephan Gruener in 2024: removed unnecessary bits for final analysis

library("gdata")
library("Seurat")
library("tidyverse")
library("Matrix")
library("ComplexHeatmap")
library("RColorBrewer")
library("clusterProfiler")
library("dorothea")
library("gplots")
library("circlize")
library("DElegate")
library("hdf5r")
library("msigdbr")

#mouse->human orthologs map (keep only mouse gene with a single unanbigous human orthologs)
m2h <- read_tsv("../data/homology_map/20200307_ensembl/mouse.txt", col_names = c("mouse", "human"))
cs <- table(m2h$mouse)
genes_single <- names(cs)[cs == 1]
m2h_single_match <- m2h %>% filter(mouse %in% genes_single)

#blacklist (genes)
bl_genes <- read_tsv("bl.txt", col_names = FALSE) %>% pull(1) %>% sort()

# prep hallmarks gene sets
hallmark <- msigdbr("Homo sapiens", "H")

# function for computing hallmark gene set module scores in a seurat object
hallmark_module_score <- function(hallmark, so, ctrl_genes, name_prefix = "") {
  for(gs in unique(hallmark$gs_name)) {
    
    gs_genes <- hallmark %>% dplyr::filter(gs_name == gs) %>% pull(gene_symbol) %>% list()
    
    so <- AddModuleScore(
      object = so,
      features = gs_genes,
      ctrl = ctrl_genes,
      name = paste0(name_prefix, gs))
      
    gc()
  }
  return(so)
}

# define PA signatures
pa_up_genes <- read_tsv(paste0(working_dir, "PA_UP.txt"), col_names = F) %>%
  pull(X1)
pa_down_genes <- read_tsv(paste0(vsc_data_dir, "PA_DOWN.txt"), col_names = F) %>%
  pullworking_dir

#PA without RPL/RPS
w <- ! grepl("^RPL|^RPS", pa_up_genes)
pa_up_genes_no_ribo <- pa_up_genes[w]
w <- ! grepl("^RPL|^RPS", pa_down_genes)
pa_down_genes_no_ribo <- pa_down_genes[w]

pa_signatures <- list("pa_up" = pa_up_genes, "pa_down" = pa_down_genes,
                      "pa_up_noRibo" = pa_up_genes_no_ribo, "pa_down_noRibo" = pa_down_genes_no_ribo)

# function for computing PA signature module scores in a seurat object
pa_module_score <- function(pa_signatures, so, ctrl_genes, name_prefix = "") {
  for(gs in names(pa_signatures)) {
    
    genes <- list(pa_signatures[[gs]])
    
    so <- AddModuleScore(
      object = so,
      features = genes,
      ctrl = ctrl_genes,
      name = paste0(name_prefix, gs))
    
    gc()
  }
  return(so)
}
