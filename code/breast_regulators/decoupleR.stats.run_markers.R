
out_folder <- "decoupleR.stats.run_markers/"
cmd <- paste0("mkdir ", out_folder)
system(cmd)


#~~~~~~~~~~~
#Ifn + Stats
############

#DimPlot(so_subset)

#FeaturePlot(so_subset, features = c("STAT1"), min.cutoff = -2, max.cutoff = 2) + 
#  scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red')

box_plt <- list()
vln_plt <- list()

#hallmarks
box_plt[["hallmarks"]] <- list()
vln_plt[["hallmarks"]] <- list()
ifn_hm <- c("post_integration_HALLMARK_INTERFERON_GAMMA_RESPONSE", 
            "post_integration_HALLMARK_INTERFERON_ALPHA_RESPONSE")
for (i in 1:length(ifn_hm)) {
  id <- ifn_hm[i]
  box_plt[["hallmarks"]][[id]] <- meta %>% 
    ggplot(aes(x = cell_type_2, y = get(id))) + 
    geom_boxplot() + 
    theme_bw() +
    ggtitle(gsub("post_integration_HALLMARK_", "", id)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Score")
  vln_plt[["hallmarks"]][[id]] <- meta %>% 
    ggplot(aes(x = cell_type_2, y = get(id))) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
    theme_bw() +
    ggtitle(gsub("post_integration_HALLMARK_", "", id)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Score") +
    ylim(c(-1, 4))
}

#signaling
box_plt[["signaling"]] <- list()
vln_plt[["signaling"]] <- list()
ifn_sig <- c("JAK-STAT")
for (i in 1:length(ifn_sig)) {
  id <- ifn_sig[i]
  box_plt[["signaling"]][[id]] <- sig_all %>% 
    dplyr::filter(source == id) %>%
    ggplot(aes(x = cell_type, y = score)) + 
    geom_boxplot() + 
    theme_bw() +
    ggtitle(id) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  vln_plt[["signaling"]][[id]] <- sig_all %>% 
    dplyr::filter(source == id) %>%
    ggplot(aes(x = cell_type, y = score)) + 
    geom_violin(trim=FALSE) +
    geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
    theme_bw() +
    ggtitle(id) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylim(c(-10,20))
}

out_file <- paste0(out_folder, "/", "boxplots.hallmarks.pdf")
pdf(out_file, width = 3.5, height = 3)
for (i in 1:length(ifn_hm)) {
  id <- ifn_hm[i]
  plot(box_plt[["hallmarks"]][[id]])
}
dev.off()

out_file <- paste0(out_folder, "/", "boxplots.signaling.pdf")
pdf(out_file, width = 3.5, height = 3)
for (i in 1:length(ifn_sig)) {
  id <- ifn_sig[i]
  plot(box_plt[["signaling"]][[id]])
}
dev.off()

out_file <- paste0(out_folder, "/", "vioplots.hallmarks.pdf")
pdf(out_file, width = 3.5, height = 3)
for (i in 1:length(ifn_hm)) {
  id <- ifn_hm[i]
  plot(vln_plt[["hallmarks"]][[id]])
}
dev.off()

out_file <- paste0(out_folder, "/", "vioplots.signaling.pdf")
pdf(out_file, width = 3.5, height = 3)
for (i in 1:length(ifn_sig)) {
  id <- ifn_sig[i]
  plot(vln_plt[["signaling"]][[id]])
}
dev.off()

