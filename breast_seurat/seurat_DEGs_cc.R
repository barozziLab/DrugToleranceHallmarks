# read seurat object
so <- readRDS(file = paste(outFolder, "so.rds", sep = ""))

#~~~~~~~~~~~~~~~~~~~
# cell cycle markers
####################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

so <- CellCycleScoring(so, 
                       s.features = s.genes, 
                       g2m.features = g2m.genes,
                       set.ident = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~
# gene sets - module scores
###########################

so <- hallmark_module_score(hallmark = hallmark, 
                            so = so, 
                            ctrl_genes = 100, 
                            name_prefix = "prior_integration_")

so <- pa_module_score(pa_signatures = pa_signatures,
                      so = so,
                      ctrl_genes = 100,
                      name_prefix = "prior_integration_")

#~~~~~~~~~~~~~
# DEG analyses
##############

print("Running DEG analysis")

options(future.globals.maxSize = 16 * 10^9)
future::plan(strategy = 'future::multicore', workers = n_cores)

## comparisons between conditions
set.seed(123)

# deseq2 method
res_deseq <- list()
for (n in 1:nrow(metadata_cmp)) {
  comp <- metadata_cmp[n,]
  res_deseq[[paste(comp$treated, "vs", comp$untreated, sep = "_")]] <-
    findDE(so, group_column = "cell_type", method = "deseq", compare = c(comp$treated, comp$untreated))
}
saveRDS(res_deseq, file = paste(outFolder, "/deg_results_deseq.rds", sep = ""))

# edgeR method
res_edger <- list()
for (n in 1:nrow(metadata_cmp)) {
  comp <- metadata_cmp[n,]
  res_edger[[paste(comp$treated, "vs", comp$untreated, sep = "_")]] <-
    findDE(so, group_column = "cell_type", method = "edger", compare = c(comp$treated, comp$untreated))
}
saveRDS(res_edger, file = paste(outFolder, "/deg_results_edger.rds", sep = ""))

#~~~~~~~~~~~~~~~~
# mean expression
#################

# save mean expression of genes across all cells to allow filtering DEGs by that
avg_vec <- apply(so@assays$RNA@data, 1, mean)

avg_expr <- tibble(gene = names(avg_vec), avg_expr = avg_vec)
write_tsv(x = avg_expr, file = paste(outFolder, "/avg_expr_breast.tsv", sep = ""))

#~~~~
#Save
#####

# save metadata
metadata_so <- so@meta.data %>% rownames_to_column(var = "rowname")
write_tsv(x = metadata_so, file = paste(outFolder, "/so_metadata_prior_integration.tsv", sep = ""))

saveRDS(so, file = paste(outFolder, "/so_degs_cc.rds", sep = ""))