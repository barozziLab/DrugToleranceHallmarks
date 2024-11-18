
library("tidyverse")
library("DESeq2")

d_complete <- read_tsv("1-s2.0-S1535610820306097-mmc3.txt")
#exclude genes from bl
bl <- read_tsv("/home/iros/imperial/datasets/UCSC/hg38_refGene_20181213.bl.txt", col_names = FALSE)
bl <- bl %>% pull(1)
d_complete <- d_complete %>% filter(! symbol %in% bl)
#average rows for the same gene
d <- d_complete %>% group_by(symbol) %>% summarize_all(mean)

d_meta <- read_tsv("1-s2.0-S1535610820306097-mmc3.meta.txt", col_types="cff")
d_meta <- d_meta %>% mutate(Sample_Treatment = paste(Sample, Treatment, sep = "-"))

d_cfr <- read_tsv("1-s2.0-S1535610820306097-mmc3.cfr.txt")

deseq_results <- list()

for (i in 1:nrow(d_cfr)) {

	id1 <- d_cfr %>% dplyr::slice(i) %>% pull("Sample_Treatment_1")
	id2 <- d_cfr %>% dplyr::slice(i) %>% pull("Sample_Treatment_2")

	coldata <- d_meta %>% filter(Sample_Treatment == id1 | Sample_Treatment == id2)
	ids <- coldata %>% pull(ID)
	counts <- d %>% select(ids) %>% as.data.frame()
	rownames(counts) <- d %>% pull(symbol)
	#round
	counts <- round(counts)
	#exclude genes with low counts
	keep <- apply(counts > 1, 1, sum) >= 2
	counts <- counts[keep,]

	dds <- DESeqDataSetFromMatrix(countData = counts,
	                              colData = coldata,
	                              design= ~ Treatment)

	colData_tmp <- colData(dds)
	colData_tmp$Treatment <- relevel(colData_tmp$Treatment, "vehicle")
	colData(dds) <- colData_tmp

	dds <- DESeq(dds)
	res <- results(dds)
	res <- res %>% as_tibble(rownames="Symbol") %>% filter(padj <= 0.05) %>% arrange(-log2FoldChange)

	id <- paste(id1, "_VS_", id2, sep = "")
	deseq_results[[id]] <- res

}

saveRDS(deseq_results, file = "deseq_results.rds")

