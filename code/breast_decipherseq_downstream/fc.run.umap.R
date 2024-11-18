
umap_nn_grid <- c(10, 15, 20, 50)
for (n_neighbors in umap_nn_grid) {
  print(paste("Running umap for", n_neighbors, "neighbors", sep = " ")) # progress statement
  gc()
  custom.config <- umap.defaults
  custom.config$n_neighbors <- n_neighbors
  umap_res <- umap_on_nmf(H.data, custom.config = custom.config)
  
  saveRDS(object = umap_res,
          file = paste0("umap_res_", n_neighbors, "_neighbors.rds"))
}

description_umap <- list()
description_split_umap <- list()
for (n_neighbors in umap_nn_grid) {
  umap_res <- readRDS(file = paste0("umap_res_", n_neighbors, "_neighbors.rds"))
  
  name <- paste("neighbors", n_neighbors, sep = "_")
  
  col <- umap_res$embeddings.labeled$color
  names(col) <- umap_res$embeddings.labeled$description
  description_umap[[name]] <- ggplot(umap_res$embeddings.labeled, aes(x = UMAP1, y = UMAP2, label = description)) +
    geom_point(aes(color = description), size = 0.1, alpha = 0.5) + 
    scale_color_manual(values = col) +
    theme_classic() +
    xlab("UMAP1") +
    ylab("UMAP2") +
    guides(colour = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
    ggtitle(name)
  
  description_split_umap[[name]] <- ggplot(umap_res$embeddings.labeled, aes(x = UMAP1, y = UMAP2, label = description)) +
    geom_point(aes(color = description), size = 0.1, alpha = 0.5) + 
    scale_color_manual(values = col) +
    theme_bw() +
    xlab("UMAP1") +
    ylab("UMAP2") +
    guides(color = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
    facet_wrap(~description) +
    ggtitle(name)
  
  ggsave(filename = paste("umap_description_res_", n_neighbors, "_neighbors.png"), description_umap[[name]], width = 18, height = 12)
  ggsave(filename = paste("umap_description_split_res_", n_neighbors, "_neighbors.png"), description_split_umap[[name]], width = 18, height = 12)
}

# create one split version to investigate how the kit version (remaining batch effects) behave

# first line, resistant, and second line treatments of each treatment line
n_neighbors <- 50
umap_res <- readRDS(file = paste0("umap_res_", n_neighbors, "_neighbors.rds"))

Lted_Tam_group <- c("MCF7_RM", "MCF7_Tam-2d", "MCF7_Lted", "MCF7_Lted_Tam-2d", "MCF7_Lted_Tam-7d")
Lted_Fulv_group <- c("MCF7_RM", "MCF7_Fulv-2d", "MCF7_Lted", "MCF7_Lted_Fulv-2d", "MCF7_Lted_Fulv-7d")

H.data_subset <- H.data %>% filter(description %in% Lted_Tam_group)
custom.config <- umap.defaults
custom.config$n_neighbors <- n_neighbors
umap_res_Tam_subset <- umap_on_nmf(H.data_subset, custom.config = custom.config)

H.data_subset <- H.data %>% filter(description %in% Lted_Fulv_group)
custom.config <- umap.defaults
custom.config$n_neighbors <- n_neighbors
umap_res_Fulv_subset <- umap_on_nmf(H.data_subset, custom.config = custom.config)


col <- umap_res_Tam_subset$embeddings.labeled$color
names(col) <- umap_res_Tam_subset$embeddings.labeled$description
description_umap[["Lted_Tam"]] <- ggplot(umap_res_Tam_subset$embeddings.labeled, aes(x = UMAP1, y = UMAP2, label = description)) +
  geom_point(aes(color = description), size = 0.1, alpha = 0.5) + 
  scale_color_manual(values = col) +
  theme_classic() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
  ggtitle("Lted_Tam")

description_split_umap[["Lted_Tam"]] <- ggplot(umap_res_Tam_subset$embeddings.labeled, aes(x = UMAP1, y = UMAP2, label = description)) +
  geom_point(aes(color = description), size = 0.1, alpha = 0.5) + 
  scale_color_manual(values = col) +
  theme_bw() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  guides(color = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
  facet_wrap(~description) +
  ggtitle("Lted_Tam")


col <- umap_res_Fulv_subset$embeddings.labeled$color
names(col) <- umap_res_Fulv_subset$embeddings.labeled$description
description_umap[["Lted_Fulv"]] <- ggplot(umap_res_Fulv_subset$embeddings.labeled, aes(x = UMAP1, y = UMAP2, label = description)) +
  geom_point(aes(color = description), size = 0.1, alpha = 0.5) + 
  scale_color_manual(values = col) +
  theme_classic() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
  ggtitle("Lted_Fulv")

description_split_umap[["Lted_Fulv"]] <- ggplot(umap_res_Fulv_subset$embeddings.labeled, aes(x = UMAP1, y = UMAP2, label = description)) +
  geom_point(aes(color = description), size = 0.1, alpha = 0.5) + 
  scale_color_manual(values = col) +
  theme_bw() +
  xlab("UMAP1") +
  ylab("UMAP2") +
  guides(color = guide_legend(override.aes = list(size = 3, shape = 16, alpha = 1))) +
  facet_wrap(~description) +
  ggtitle("Lted_Fulv")

