#Load Seurat object
inFile <- paste(outFolder, "/so_degs_cc.rds", sep = "")
so.cca <- readRDS(inFile)

#~~~
#CCA
####

#need to increase this parameter from 500M to 125000M (125000*1024^2)
options(future.globals.maxSize = 131072000000)

seurat_obj_cca.list.complete <- SplitObject(so.cca, split.by = "kit")

for (i in 1:length(seurat_obj_cca.list.complete)) {
  
  seurat_obj_cca.list.complete[[i]] <- SCTransform(seurat_obj_cca.list.complete[[i]], verbose = FALSE, variable.features.n = n_features)
  
  #qc UMAP
  
  seurat_obj_cca.tmp <- RunPCA(seurat_obj_cca.list.complete[[i]],
                               verbose = FALSE,
                               npcs = 100)
  seurat_obj_cca.tmp <- RunUMAP(seurat_obj_cca.tmp, dims = 1:50)
  
  p <- DimPlot(object = seurat_obj_cca.tmp,
               reduction = 'umap',
               group.by = 'cell_type_id',
               label = TRUE,
               repel = TRUE,
               label.size = 2)
  
  if (names(seurat_obj_cca.list.complete)[i] == "v2") {
    width <- 8
  } else if (names(seurat_obj_cca.list.complete)[i] == "v3") {
    width <- 12
  }
  
  outFile <- paste(outFolder, "/data.umap.", names(seurat_obj_cca.list.complete)[i], ".pdf", sep = "")
  pdf(outFile, width = width, height = 5)
  plot(p)
  dev.off()
  
}

seurat_obj_cca.list <- seurat_obj_cca.list.complete

seurat_obj_cca.features <- SelectIntegrationFeatures(object.list = seurat_obj_cca.list, nfeatures = n_features)

seurat_obj_cca.list <- PrepSCTIntegration(object.list = seurat_obj_cca.list,
                                          anchor.features = seurat_obj_cca.features,
                                          verbose = FALSE)

seurat_obj_cca.anchors <- FindIntegrationAnchors(object.list = seurat_obj_cca.list,
                                                 normalization.method = "SCT",
                                                 anchor.features = seurat_obj_cca.features,
                                                 verbose = FALSE)

seurat_obj_cca.integrated <- IntegrateData(anchorset = seurat_obj_cca.anchors,
                                           normalization.method = "SCT",
                                           verbose = FALSE)

seurat_obj_cca.integrated <- RunPCA(seurat_obj_cca.integrated,
                                    verbose = FALSE,
                                    npcs = 100)
#plot(seurat_obj_cca.integrated@reductions$pca@stdev)
pca_dim_sel <- 50

seurat_obj_cca.integrated <- RunUMAP(seurat_obj_cca.integrated, dims = 1:pca_dim_sel)

seurat_obj_cca.integrated <- FindNeighbors(seurat_obj_cca.integrated, dims = 1:pca_dim_sel)

# clustering using leiden algorithm
seurat_obj_cca.integrated <- FindClusters(seurat_obj_cca.integrated, resolution = clustering_res, algorithm = 4, method = "igraph")

p <- DimPlot(object = seurat_obj_cca.integrated,
             reduction = 'umap',
             group.by = 'cell_type_id',
             label = TRUE,
             repel = TRUE,
             label.size = 2)

outFile <- paste(outFolder, "/UMAP.cell_type_id.pdf", sep = "")
pdf(outFile, width = 16, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = seurat_obj_cca.integrated,
             reduction = 'umap',
             group.by = 'cell_type',
             label = TRUE,
             repel = TRUE,
             label.size = 2)

outFile <- paste(outFolder, "/UMAP.cell_type.pdf", sep = "")
pdf(outFile, width = 10, height = 5)
plot(p)
dev.off()

p <- DimPlot(object = seurat_obj_cca.integrated,
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

p <- DimPlot(object = seurat_obj_cca.integrated,
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

seurat_obj_cca.integrated <- CellCycleScoring(seurat_obj_cca.integrated, 
                       s.features = s.genes, 
                       g2m.features = g2m.genes,
                       set.ident = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~
# gene sets - module scores
###########################

seurat_obj_cca.integrated <- hallmark_module_score(hallmark = hallmark, 
                                                   so = seurat_obj_cca.integrated, 
                                                   ctrl_genes = 100, 
                                                   name_prefix = "post_integration_")

seurat_obj_cca.integrated <- pa_module_score(pa_signatures = pa_signatures,
                                             so = seurat_obj_cca.integrated,
                                             ctrl_genes = 100,
                                             name_prefix = "post_integration_")

#~~~~
#Save
#####

# save metadata
metadata_so_int <- seurat_obj_cca.integrated@meta.data %>% rownames_to_column(var = "rowname")
write_tsv(x = metadata_so_int, file = paste(outFolder, "/so_metadata_post_integration.tsv", sep = ""))

# save integrated seurat object
saveRDS(seurat_obj_cca.integrated, file = paste(outFolder, "/so_cca_integrated_cc.rds", sep = ""))
