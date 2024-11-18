
library("tidyverse")
library("GEOquery")
library("lumi")
library("limma")

gset <- getGEO("GSE102124", GSEMatrix = TRUE, getGPL = TRUE)

expr <- exprs(gset$GSE102124_series_matrix.txt.gz)

meta_genes <- fData(gset$GSE102124_series_matrix.txt.gz)
meta_genes_tib <- meta_genes %>% as_tibble()

#subset on genes
w <- ! is.na(meta_genes_tib$ID)

meta_genes_tib <- meta_genes_tib[w,]
expr <- expr[w,]

#extract gene symbols

meta_genes_tib <- meta_genes_tib %>% 
  rowwise() %>% 
  mutate(gene_symbol = strsplit(gene_assignment, "//")[[1]][2]) %>%
  ungroup() %>%
  mutate(gene_symbol = gsub(" ", "", gene_symbol))

#exclude --- and OTT/LOC/

w <- ! grepl("^LOC|^OTT|^---", meta_genes_tib$gene_symbol)

meta_genes_tib <- meta_genes_tib[w,]
expr <- expr[w,]

meta_samples <- pData(gset$GSE102124_series_matrix.txt.gz)

#normalize with lumi [not needed, see below]
#expr_norm <- smoothQuantileNormalization(expr, ref = expr[,"GSM2724470"], logMode=FALSE)

#expr already log2

#boxplot(expr)
#it looks like they already are normalized

#exclude sd = 0
sd_gene <- apply(expr, 1, sd)
expr <- expr[sd_gene > 0,]
meta_genes_tib <- meta_genes_tib[sd_gene > 0,]

#DEGs using limma

treat <- meta_samples[,"treatment group:ch1"]

design <- model.matrix(~0+treat)
fit <- lmFit(expr, design)
contrasts <- makeContrasts(treattreated-treatuntreated, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)

#
res <- topTable(fit2, coef=1, p.value=0.05, number=20000) %>% as_tibble(rownames="ID")
res <- res %>% 
  left_join(meta_genes_tib %>% select(ID, gene_symbol), by = "ID") %>% 
  select(ID, gene_symbol, logFC, AveExpr, P.Value, adj.P.Val)

write_tsv(res, file = "limma_results.txt")

