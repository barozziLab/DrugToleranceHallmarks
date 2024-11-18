
source("downstream.libs.R")

outFolder <- "downstream.pervasiveness.results/"
cmd <- paste("mkdir", outFolder)
system(cmd)

cc.genes <- readRDS("../data/cc.genes_seurat.rds")

cc.genes.vec <- cc.genes %>% unlist() %>% as.vector()

#Hallmark gene sets
m_cat <- msigdbr(species = "Homo sapiens", category = "H")
m_t2g <- m_cat %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()


#~~~~~~~~~~~~~~
#Neo Signatures
###############

#published
source("downstream.published.dataPrep.R")

#ours
sleuth_table <- read_tsv("sleuth.results/sleuth.results.txt")
gs_neo[["HR_novel_UP"]] <- sleuth_table %>% dplyr::filter(qval <= 0.05 & b >= 0.5) %>% pull(target_id) %>% sort()
gs_neo[["HR_novel_DOWN"]] <- sleuth_table %>% dplyr::filter(qval <= 0.05 & b <= -0.5) %>% pull(target_id) %>% sort()

#convert to t2g
neo_t2g <- c()
for (gs_name in names(gs_neo)) {
  neo_t2g <- rbind(neo_t2g, cbind(gs_name, gs_neo[[gs_name]]))
}
neo_t2g <- data.frame(gs_name = neo_t2g[,1], gene_symbol = neo_t2g[,2])

#Exclude Ribosomal genes
w <- ! grepl("^RPL|^RPS", neo_t2g$gene_symbol)
neo_t2g <- neo_t2g[w,]

#Exclude CC genes
w <- ! neo_t2g$gene_symbol %in% cc.genes.vec
neo_t2g <- neo_t2g[w,]


#~~~~~~~~~~~~~~
#DTP Signatures
###############

dtp_ds <- read_tsv("../literature_meta-analysis/DTP/DTPs.markers.txt")

#Exclude Diapause & CTC
dtp_ds <- dtp_ds %>% 
  dplyr::filter(!grepl("Diapause|CTC|common|both", gs_name))

dtp_t2g <- dtp_ds %>% 
  dplyr::filter(species == "human") %>% 
  dplyr::select(gs_name, gene) %>% 
  dplyr::rename(gene_symbol = gene) %>% 
  as.data.frame() %>%
  unique()

#Exclude Ribosomal genes
w <- ! grepl("^RPL|^RPS", dtp_t2g$gene_symbol)
dtp_t2g <- dtp_t2g[w,]

#ONLY APPLIES TO MOUSE
#Exclude Ceacam/Psg genes (some issues/inflation with paralogs mapping)
#w <- ! grepl("^Ceacam|^Psg", dtp_t2g$gene_symbol)
#dtp_t2g <- dtp_t2g[w,]

#Exclude CC genes
w <- ! dtp_t2g$gene_symbol %in% cc.genes.vec
dtp_t2g <- dtp_t2g[w,]


#~~~~~~~~~~~~~~~~~~~
#Pervasiveness - Neo
####################

#Unique signature considered
neo_perv_sig_id <- gsub("_UP|_DOWN", "", neo_t2g[grepl("_UP|_DOWN", neo_t2g$gs_name), "gs_name"]) %>% unique()
out_F <- paste(outFolder, "/NEO-signatures.ids.txt", sep = "")
write_tsv(x = neo_perv_sig_id %>% as_tibble(), file = out_F, col_names = FALSE)
#and corresponding markers
neo_perv_sig_markers <- neo_t2g %>% dplyr::filter(gs_name %in% neo_t2g[grepl("_UP|_DOWN", neo_t2g$gs_name), "gs_name"])
out_F <- paste(outFolder, "/NEO-signatures.markers.txt", sep = "")
write_tsv(x = neo_perv_sig_markers %>% as_tibble(), file = out_F)

#Derive up-regulated genes shared in 2+ (published) signatures
neo_perv_up <- neo_t2g[grepl("_UP", neo_t2g$gs_name),] %>% unique() %>% group_by(gene_symbol) %>% summarise(n = n()) %>% arrange(-n)

#Same for down-regulated
neo_perv_down <- neo_t2g[grepl("_DOWN", neo_t2g$gs_name),] %>% unique() %>% group_by(gene_symbol) %>% summarise(n = n()) %>% arrange(-n)

#Merge
neo_perv <- rbind(neo_perv_up %>% mutate(class = "up"), neo_perv_down %>% mutate(class = "down")) %>% arrange(-n)

#Merge
neo_perv <- rbind(neo_perv_up %>% mutate(class = "up"), neo_perv_down %>% mutate(class = "down")) %>% arrange(-n)
neo_perv <- neo_perv %>% 
  mutate(class_2 = ifelse(n >= 4 & class == "up", "highly_pervasive_up", ifelse(n >= 4 & class == "down", "highly_pervasive_down", "not_pervasive")))
out_F <- paste(outFolder, "/NEO-signatures.txt", sep = "")
write_tsv(x = neo_perv %>% as_tibble(), file = out_F)

#GSEA on top pervasive

i <- 0
g_names <- c()
genes_ns <- c()
res_all <- c()

for (cl in c("highly_pervasive_up", "highly_pervasive_down")) {
  
  i <- i + 1
  
  genes <- neo_perv %>% filter(class_2 == cl) %>% pull(gene_symbol) %>% sort()
  
  if (length(genes) >= 10) {
    
    res <- enricher(gene = genes, TERM2GENE = m_t2g, maxGSSize = 1500)
    
    if (nrow(res) > 0) {
      
      g_names <- c(g_names, cl)
      genes_ns <- c(genes_ns, length(genes))
      
      #keep only those significant
      keep <- res$p.adjust <= 0.05
      res_keep <- res[keep,]
      
      #append to master table
      res_keep_mt <- res_keep
      res_keep_mt$cluster_id <- cl
      
      res_all <- rbind(res_all, res_keep_mt)
      
      if (i == 1) {
        m_cat_res_padj <- res_keep[,c("ID", "p.adjust")]
        m_cat_res_n <- res_keep[,c("ID", "Count")]
      } else {
        m_cat_res_padj <- merge(m_cat_res_padj, res_keep[,c("ID", "p.adjust")], by="ID", all.x=TRUE, all.y=TRUE)
        m_cat_res_n <- merge(m_cat_res_n, res_keep[,c("ID", "Count")], by="ID", all.x=TRUE, all.y=TRUE)
      }
      
    }
    
  }
  
}

rownames(m_cat_res_padj) <- m_cat_res_padj[,"ID"]
m_cat_res_padj <- m_cat_res_padj[,-1]
colnames(m_cat_res_padj) <- g_names
m_cat_res_padj[is.na(m_cat_res_padj)] <- 1

rownames(m_cat_res_n) <- m_cat_res_n[,"ID"]
m_cat_res_n <- m_cat_res_n[,-1]
colnames(m_cat_res_n) <- g_names
m_cat_res_n[is.na(m_cat_res_n)] <- 0

m_cat_res_frac <- t(t(m_cat_res_n) / genes_ns)

#Save full results

out_F <- paste(outFolder, "/NEO-signatures.gsea.txt", sep = "")
write_tsv(x = res_all %>% as_tibble(), file = out_F)

#Bubble-plot

frac <- m_cat_res_frac %>% as.data.frame() %>% rownames_to_column(var = "id") %>% as_tibble()
padj <- m_cat_res_padj %>% as.data.frame() %>% rownames_to_column(var = "id") %>% as_tibble()

frac_long <- frac %>% pivot_longer(cols = -c("id"), names_to = "group", values_to = "fraction")
padj_long <- padj %>% pivot_longer(cols = -c("id"), names_to = "group", values_to = "padjust")

data <- frac_long %>% left_join(padj_long, by = c("id", "group"))

padj_cl <- rep("ns", nrow(data))
padj_cl[data$padjust <= 0.05] <- "adjusted_pval <= 0.05"
padj_cl[data$padjust <= 1e-5] <- "adjusted_pval <= 1e-5"
padj_cl[data$padjust <= 1e-10] <- "adjusted_pval <= 1e-10"
padj_cl <- factor(padj_cl, levels = rev(c("ns", "adjusted_pval <= 0.05", "adjusted_pval <= 1e-5", "adjusted_pval <= 1e-10")))

data <- data %>% mutate(padjust_cl = padj_cl)

#re-order groups
data$group <- factor(data$group, levels = g_names)

#re-order rows according to delta (up-down fraction)
w_ord <- order(m_cat_res_frac[,"highly_pervasive_up"] - m_cat_res_frac[,"highly_pervasive_down"])
data$id <- factor(data$id, levels = rownames(m_cat_res_frac)[w_ord])

plt <- ggplot(data, aes(x = group, y = id)) + 
  geom_point(aes(size = fraction, fill = padjust_cl), alpha = 0.75, shape = 21) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_brewer(palette = "Reds", direction = -1) +
  scale_size(range = c(.1, 5), name="Fraction of Genes")

out_F <- paste(outFolder, "/NEO-signatures.gsea.bubbleplot.pdf", sep = "")
pdf(out_F, width = 6.5, height = 6)
plot(plt)
dev.off()


#~~~~~~~~~~~~~~~~~~~
#Pervasiveness - DTP
####################

#Unique signature considered
dtp_perv_sig_id <- gsub("-up|-down", "", dtp_t2g[grepl("-up|-down", dtp_t2g$gs_name), "gs_name"]) %>% unique()
out_F <- paste(outFolder, "/DTP-signatures.ids.txt", sep = "")
write_tsv(x = dtp_perv_sig_id %>% as_tibble(), file = out_F, col_names = FALSE)
#and corresponding markers
dtp_perv_sig_markers <- dtp_t2g %>% dplyr::filter(gs_name %in% dtp_t2g[grepl("-up|-down", dtp_t2g$gs_name), "gs_name"])
out_F <- paste(outFolder, "/DTP-signatures.markers.txt", sep = "")
write_tsv(x = dtp_perv_sig_markers %>% as_tibble(), file = out_F)

#Derive up-regulated genes shared in 2+ (published) signatures
dtp_perv_up <- dtp_t2g[grepl("-up", dtp_t2g$gs_name),] %>% unique() %>% group_by(gene_symbol) %>% summarise(n = n()) %>% arrange(-n)

#Same for down-regulated
dtp_perv_down <- dtp_t2g[grepl("-down", dtp_t2g$gs_name),] %>% unique() %>% group_by(gene_symbol) %>% summarise(n = n()) %>% arrange(-n)

#Merge
dtp_perv <- rbind(dtp_perv_up %>% mutate(class = "up"), dtp_perv_down %>% mutate(class = "down")) %>% arrange(-n)
dtp_perv <- dtp_perv %>% 
  mutate(class_2 = ifelse(n >= 4 & class == "up", "highly_pervasive_up", ifelse(n >= 4 & class == "down", "highly_pervasive_down", "not_pervasive")))
out_F <- paste(outFolder, "/DTP-signatures.txt", sep = "")
write_tsv(x = dtp_perv %>% as_tibble(), file = out_F)

#GSEA on top pervasive

i <- 0
g_names <- c()
genes_ns <- c()
res_all <- c()

for (cl in c("highly_pervasive_up", "highly_pervasive_down")) {
  
  i <- i + 1
  
  genes <- dtp_perv %>% filter(class_2 == cl) %>% pull(gene_symbol) %>% sort()
  
  if (length(genes) >= 10) {
    
    res <- enricher(gene = genes, TERM2GENE = m_t2g, maxGSSize = 1500)
    
    if (nrow(res) > 0) {
      
      g_names <- c(g_names, cl)
      genes_ns <- c(genes_ns, length(genes))
      
      #keep only those significant
      keep <- res$p.adjust <= 0.05
      res_keep <- res[keep,]
      
      #append to master table
      res_keep_mt <- res_keep
      res_keep_mt$cluster_id <- cl
      
      res_all <- rbind(res_all, res_keep_mt)
      
      if (i == 1) {
        m_cat_res_padj <- res_keep[,c("ID", "p.adjust")]
        m_cat_res_n <- res_keep[,c("ID", "Count")]
      } else {
        m_cat_res_padj <- merge(m_cat_res_padj, res_keep[,c("ID", "p.adjust")], by="ID", all.x=TRUE, all.y=TRUE)
        m_cat_res_n <- merge(m_cat_res_n, res_keep[,c("ID", "Count")], by="ID", all.x=TRUE, all.y=TRUE)
      }
      
    }
    
  }
  
}

rownames(m_cat_res_padj) <- m_cat_res_padj[,"ID"]
m_cat_res_padj <- m_cat_res_padj[,-1]
colnames(m_cat_res_padj) <- g_names
m_cat_res_padj[is.na(m_cat_res_padj)] <- 1

rownames(m_cat_res_n) <- m_cat_res_n[,"ID"]
m_cat_res_n <- m_cat_res_n[,-1]
colnames(m_cat_res_n) <- g_names
m_cat_res_n[is.na(m_cat_res_n)] <- 0

m_cat_res_frac <- t(t(m_cat_res_n) / genes_ns)

#Save full results

out_F <- paste(outFolder, "/DTP-signatures.gsea.txt", sep = "")
write_tsv(x = res_all %>% as_tibble(), file = out_F)

#Bubble-plot

frac <- m_cat_res_frac %>% as.data.frame() %>% rownames_to_column(var = "id") %>% as_tibble()
padj <- m_cat_res_padj %>% as.data.frame() %>% rownames_to_column(var = "id") %>% as_tibble()

frac_long <- frac %>% pivot_longer(cols = -c("id"), names_to = "group", values_to = "fraction")
padj_long <- padj %>% pivot_longer(cols = -c("id"), names_to = "group", values_to = "padjust")

data <- frac_long %>% left_join(padj_long, by = c("id", "group"))

padj_cl <- rep("ns", nrow(data))
padj_cl[data$padjust <= 0.05] <- "adjusted_pval <= 0.05"
padj_cl[data$padjust <= 1e-5] <- "adjusted_pval <= 1e-5"
padj_cl[data$padjust <= 1e-10] <- "adjusted_pval <= 1e-10"
padj_cl <- factor(padj_cl, levels = rev(c("ns", "adjusted_pval <= 0.05", "adjusted_pval <= 1e-5", "adjusted_pval <= 1e-10")))

data <- data %>% mutate(padjust_cl = padj_cl)

#re-order groups
data$group <- factor(data$group, levels = g_names)

#re-order rows according to delta (up-down fraction)
w_ord <- order(m_cat_res_frac[,"highly_pervasive_up"] - m_cat_res_frac[,"highly_pervasive_down"])
data$id <- factor(data$id, levels = rownames(m_cat_res_frac)[w_ord])

plt <- ggplot(data, aes(x = group, y = id)) + 
  geom_point(aes(size = fraction, fill = padjust_cl), alpha = 0.75, shape = 21) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_brewer(palette = "Reds", direction = -1) +
  scale_size(range = c(.1, 5), name="Fraction of Genes")

out_F <- paste(outFolder, "/DTP-signatures.gsea.bubbleplot.pdf", sep = "")
pdf(out_F, width = 6.5, height = 7)
plot(plt)
dev.off()

##############

#~~~~~
#Plots
######

#Pervasiveness master table

perv_all <- rbind(dtp_perv_up %>% mutate(group = "DTP_up") %>% mutate(frac = n / max(n)), 
      dtp_perv_down %>% mutate(group = "DTP_down") %>% mutate(frac = n / max(n)),
      neo_perv_up %>% mutate(group = "Neo_up") %>% mutate(frac = n / max(n)),
      neo_perv_down %>% mutate(group = "Neo_down") %>% mutate(frac = n / max(n)))

plt <- ggplot(perv_all, aes(x = n)) + 
    geom_histogram() +
    facet_wrap(~group) + 
    scale_y_log10() +
    xlab("Pervasiveness") +
    ylab("") +
    scale_x_continuous(breaks = c(1,3,5,7,9)) +
    theme_bw()

out_F <- paste(outFolder, "/histograms.pdf", sep = "")
pdf(out_F, width = 5, height = 3)
plot(plt)
dev.off()

