
m_t2g <- list()

#Dorothea (only A-D confidence)
dorothea_hs_f <- dorothea_hs %>% dplyr::filter(confidence != "E")
m_t2g[["Dorothea"]] <- data.frame(gs_name = dorothea_hs_f %>% pull(tf),
                                  gene_symbol = dorothea_hs_f %>% pull(target),
                                  stringsAsFactors = FALSE)

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

#Meta-programs Tirosh
MP <- readRDS("../data/MP.rds")
m_t2g[["MP"]] <- MP$hs %>% dplyr::rename(gene_symbol = gene) %>% as.data.frame()

#Meta-programs Tirosh + PA
m_t2g[["MP_and_PA"]] <- rbind(m_t2g[["PA"]], m_t2g[["MP"]])

#neo-adjuvant ET
neoadj_data <- read_tsv("../primary_data_neo-ET/sleuth.results.txt")
neoadj_up <- neoadj_data %>% dplyr::filter(qval <= 0.05 & b >= 1) %>% pull(target_id) %>% sort()
neoadj_dn <- neoadj_data %>% dplyr::filter(qval <= 0.05 & b <= -1) %>% pull(target_id) %>% sort()
m_t2g[["Neo_ET"]] <- rbind(cbind("Neo_ET_Up", neoadj_up), cbind("Neo_ET_Down", neoadj_dn))
m_t2g[["Neo_ET"]] <- as.data.frame(m_t2g[["Neo_ET"]])
colnames(m_t2g[["Neo_ET"]]) <- c("gs_name", "gene_symbol")

#adult breast
adult_breast <- read_tsv("../data/breast_adult_kumar_et_al_2023/markers.txt")
m_t2g[["Breast_Adult"]] <- adult_breast %>% 
  arrange(gs_name, gene) %>%
  dplyr::filter(species == "hs") %>%
  dplyr::select(gs_name, gene) %>%
  dplyr::rename(gene_symbol = gene) %>%
  as.data.frame()

#filter m_t2g for the actual analysis
w <- names(m_t2g) %in% c("Hallmarks", "Pathways", "MP", "Dorothea", "Neo_ET", "Breast_Adult")
m_t2g_f <- m_t2g[w]
