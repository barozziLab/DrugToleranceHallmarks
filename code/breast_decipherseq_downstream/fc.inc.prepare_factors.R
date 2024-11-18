
#NB: this is W (!H in this case)

program_ids <- rownames(nmf_res$MCF7$W)
program_ids <- str_remove(program_ids, "R33_Program")

W <- nmf_res$MCF7$W %>% t() %>% as.data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble()
colnames(W) <- str_remove(colnames(W), "R33_Program")

top_genes_n <- 200

mrk_list <- list()
for (program_id in program_ids) {
  mrk_list[[program_id]] <- W %>% arrange(-get(program_id)) %>% head(top_genes_n) %>% pull(gene)
}

mrk_tib <- c()
for (program_id in names(mrk_list)) {
  mrk_tib <- rbind(mrk_tib, cbind(program_id, mrk_list[[program_id]] ))
}
mrk_tib <- tibble(program = mrk_tib[,1], 
                  gene_symbol = mrk_tib[,2])
