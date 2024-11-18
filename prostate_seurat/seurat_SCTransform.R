#Load Seurat object
inFile <- paste(outFolder, "/so_degs_cc.rds", sep = "")
so <- readRDS(inFile)

#~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCTransform normalization
###########################

so <- SCTransform(so, verbose = FALSE, variable.features.n = n_features)

so <- RunPCA(so,
             verbose = FALSE,
             npcs = 100)

pca_dim_sel <- 50

so <- RunUMAP(so, dims = 1:pca_dim_sel)

so <- FindNeighbors(so, dims = 1:pca_dim_sel)

# clustering using leiden algorithm
so <- FindClusters(so, resolution = clustering_res, algorithm = 4, method = "igraph")

p <- DimPlot(object = so,
             reduction = 'umap',
             group.by = 'cell_type_id',
             label = TRUE,
             repel = TRUE,
             label.size = 2)

outFile <- paste(outFolder, "/UMAP.cell_type_id.pdf", sep = "")
pdf(outFile, width = 16, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so,
             reduction = 'umap',
             group.by = 'cell_type',
             label = TRUE,
             repel = TRUE,
             label.size = 2)

outFile <- paste(outFolder, "/UMAP.cell_type.pdf", sep = "")
pdf(outFile, width = 10, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so,
             reduction = 'umap',
             group.by = 'cell_type',
             split.by = 'kit',
             ncol = 2,
             label = TRUE,
             repel = TRUE)

outFile <- paste(outFolder, "/UMAP.cell_type.split_by=kit.pdf", sep = "")
pdf(outFile, width = 15, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = so,
             reduction = 'umap',
             group.by = 'kit',
             split.by = 'cell_type',
             ncol = 4)

outFile <- paste(outFolder, "/UMAP.kit.split_by=cell_type.pdf", sep = "")
pdf(outFile, width = 10, height = 20)
plot(p)
dev.off()

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
                            name_prefix = "post_normalization_")

so <- pa_module_score(pa_signatures = pa_signatures,
                      so = so,
                      ctrl_genes = 100,
                      name_prefix = "post_normalization_")

#~~~~
#Save
#####

# save metadata
metadata_so_int <- so@meta.data %>% rownames_to_column(var = "rowname")
write_tsv(x = metadata_so_int, file = paste(outFolder, "/so_metadata_post_normalization.tsv", sep = ""))

# save integrated seurat object
saveRDS(so, file = paste(outFolder, "/so_normalized_cc.rds", sep = ""))
