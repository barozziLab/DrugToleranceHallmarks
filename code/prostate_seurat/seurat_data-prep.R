
#~~~~~~~~~~~~~~
#Merge Datasets
###############

for (i in 1:nrow(metadata)) {
  
  #load
  sample_name <- metadata[i, "id"]
  h5_file_path <- paste(metadata[i, "stats_folder"], metadata[i, "stats_file"], sep = "")
  d <- Read10X_h5(h5_file_path)
  
  #subsample
  if (!is.na(sampling_N)) {
    sampling_N_real <- min(sampling_N, ncol(d))
    set.seed(i)
    keep <- sample(1:ncol(d), sampling_N_real)
    d_s <- d[, keep]
  } else {
    d_s <- d
  }
  colnames(d_s) <- paste(sample_name, colnames(d_s), sep = "_")
  md_s <- data.frame(cell_type_id = rep(sample_name, ncol(d_s)),
                     batch = metadata[i, "batch"],
                     kit = metadata[i, "kit"],
                     cell_type = metadata[i, "description"],
                     stringsAsFactors = FALSE)
  rownames(md_s) <- colnames(d_s)
  
  #exclude blacklisted genes
  keep <- !rownames(d_s) %in% bl_genes
  d_s <- d_s[keep,]
  
  #create or append to seurat object
  if (i == 1) {
    so.combined <- CreateSeuratObject(counts = d_s,
                                      project = sample_name,
                                      meta.data = md_s)
  } else {
    so <- CreateSeuratObject(counts = d_s,
                             project = sample_name,
                             meta.data = md_s)
    so.tmp <- merge(so.combined, y = so, project = "ICL")
    so.combined <- so.tmp
  }
  
}

#force orig.ident to be equal to the sample
so.combined@meta.data$orig.ident <- so.combined@meta.data$cell_type_id

rm("d", "so", "so.tmp", "d_s", "md_s")

#add mito fraction
so.combined[["percent.mt"]] <- PercentageFeatureSet(so.combined, pattern = "^MT-")

#qc boxplots
out_path <- paste(outFolder, "data.qc.unfiltered.png", sep = "")
png(out_path, width = 800, height = 1800)
p <-
  VlnPlot(
    so.combined,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 1,
    pt.size = 0.1,
    group.by = "cell_type_id"
  )
plot(p)
dev.off()

###############

#~~~~~~~~~~~~~~~~~~~~~~
#Doublet Identification
#######################

#Must run SEPARATELY for each CAPTURE

#FROZEN, IT REQUIRES R>=4.0
#library(scDblFinder)

#See http://www.bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/4_computeDoubletDensity.html
#How to turn score to classification? Plot distribution?
#Add score and classification to metadata

#######################

#~~~~~~~~~
#Filtering
##########

#TODO: ADD DOUBLET FILTERING

so.combined <-
  subset(
    so.combined,
    percent.mt <= percent.mt_max &
      nFeature_RNA >= nFeature_RNA_min &
      nCount_RNA >= nCount_RNA_min &
      nCount_RNA <= nCount_RNA_max
  )

out_path <- paste(outFolder, "data.qc.filtered.png", sep = "")
png(out_path, width = 800, height = 1800)
p <-
  VlnPlot(
    so.combined,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 1,
    pt.size = 0.1,
    group.by = "cell_type_id"
  )
plot(p)
dev.off()

#exclude mito genes and set min.features to 100
keep <- grep("^MT-", rownames(so.combined@assays$RNA), invert = TRUE)
so.combined.filtered <- CreateSeuratObject(counts = so.combined@assays$RNA[keep,],
                                           project = sample_name,
                                           meta.data = so.combined@meta.data,
                                           min.features = 100)

##########

#~~~~
#Save
#####

outFile <- paste(outFolder, "/so.rds", sep = "")
saveRDS(so.combined.filtered, file = outFile)

####
