
library("tidyverse")
library("GEOquery")
library("lumi")
library("limma")

gse <- getGEO('GSE59515', GSEMatrix=TRUE)

exp <- exprs(gse$GSE59515_series_matrix.txt.gz)

meta_genes <- fData(gse$GSE59515_series_matrix.txt.gz)
meta_genes_tib <- meta_genes[,c("ID","Symbol")] %>% as_tibble()

meta_samples <- pData(gse$GSE59515_series_matrix.txt.gz)

tp <- meta_samples[,"time point:ch1"]
tp[tp == "pre-treatment"] <- "_pt"
tp[tp == "3 months of letrozole treatment"] <- "_3mo"
tp[tp == "2 wks of letrozole treatment"] <- "_2w"

patient <- meta_samples[,"patient id:ch1"]

#normalize with lumi [not needed, see below]
#exp_norm <- smoothQuantileNormalization(exp, ref = exp[,"GSM1438585"], logMode=FALSE)

#negative values to zeros, then log2
exp[exp < 1] <- 1
exp <- log2(exp)

#boxplot(exp)
#it looks like they already are normalized

#exclude sd = 0
sd_gene <- apply(exp, 1, sd)
exp <- exp[sd_gene > 0,]

#exclude lowly expressed genes [hard to tell, info about pval detection lost]

#DEGs using limma [separately 2 weeks and 3 months]
design <- model.matrix(~0+tp+patient)
fit <- lmFit(exp, design)
contrasts <- makeContrasts(tp_2w-tp_pt, tp_3mo-tp_pt, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2, trend=TRUE)

#2w
d_2w <- topTable(fit2, coef=1, p.value=0.05, number=20000) %>% as_tibble(rownames="ID")
d_2w <- d_2w %>% left_join(meta_genes_tib, by = "ID") %>% select(ID, Symbol, logFC, AveExpr, P.Value, adj.P.Val)

#3mo
d_3mo <- topTable(fit2, coef=2, p.value=0.05, number=20000) %>% as_tibble(rownames="ID")
d_3mo <- d_3mo %>% left_join(meta_genes_tib, by = "ID") %>% select(ID, Symbol, logFC, AveExpr, P.Value, adj.P.Val)

limma_results <- list("2w"=d_2w, "3mo"=d_3mo)

saveRDS(limma_results, file = "limma_results.rds")

