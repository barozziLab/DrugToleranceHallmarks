
m_t2g <- list()

#Msigdb hallmarks
m_cat <- msigdbr(species = "Homo sapiens", category = "H")
m_t2g[["Hallmarks"]] <- m_cat %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

#Msigdb pathways
m_cat <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
m_t2g[["Pathways"]] <- m_cat %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()

#PA
pa_up_genes <- read_tsv("../data/PA_UP.txt", col_names = FALSE) %>% pull(1) %>% sort()
m_t2g[["PA"]] <- cbind("PA_UP", pa_up_genes)
pa_down_genes <- read_tsv("../data/PA_DOWN.txt", col_names = FALSE) %>% pull(1) %>% sort()
m_t2g[["PA"]] <- rbind(m_t2g[["PA"]], cbind("PA_DOWN", pa_down_genes))

#PA without RPL/RPS
w <- ! grepl("^RPL|^RPS", pa_up_genes)
m_t2g[["PA"]] <- rbind(m_t2g[["PA"]], cbind("PA_UP_noRibo", pa_up_genes[w]))
w <- ! grepl("^RPL|^RPS", pa_down_genes)
m_t2g[["PA"]] <- rbind(m_t2g[["PA"]], cbind("PA_DOWN_noRibo", pa_down_genes[w]))

#Finalise PA
m_t2g[["PA"]] <- data.frame(m_t2g[["PA"]])
colnames(m_t2g[["PA"]]) <- colnames(m_t2g[["Pathways"]])

