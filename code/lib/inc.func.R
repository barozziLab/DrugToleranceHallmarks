
umap_on_nmf <- function(H.data, custom.config = NULL) {
  
  data <- H.data[, grep("^R[0-9]*_", colnames(H.data))]
  metadata <- H.data[, grep("^R[0-9]*_", colnames(H.data), invert = TRUE)]
  
  if(!is.null(custom.config)) {
    H.umap <- umap(data, config = custom.config)
  } else {
    H.umap <- umap(data)
  }
  
  embeddings.labeled <- cbind(metadata, H.umap$layout) %>%
    rename(UMAP1 = "1", UMAP2 = "2") %>%
    as_tibble()
  
  umap.plot <- ggplot(embeddings.labeled, aes(UMAP1, UMAP2, color = sample)) +
    geom_point(size = 0.2) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_classic()
  
  return(list(H.umap = H.umap, 
              embeddings.labeled = embeddings.labeled,
              umap.plot = umap.plot))
}

jaccard <- function(a, b) {
  
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  
  return (intersection/union)
  
}
