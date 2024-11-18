
source("sleuth.libs.R")

source("sleuth.gene_sets.R")

#load Sleuth results
inF <- paste(outFolder, "/sleuth.results.txt", sep="")
sleuth_table <- read_tsv(file = inF)


#~~~~~~~~~~~~~~~~~~~
#GSEA on the results
####################

# Prepare ranked list of genes
# Rank by -log10(q) * sign*(b)
ranked_genes <- sleuth_table %>% 
                  mutate(score = -log10(qval)*sign(b)) %>%
                  arrange(score)
ranked_genes_lst <- ranked_genes %>% pull(score)
names(ranked_genes_lst) <- ranked_genes %>% pull(target_id)

# GSEA
fgseaRes <- list()
for (set_id in names(gs_lst)) {
  fgseaRes[[set_id]] <- fgsea(pathways = gs_lst[[set_id]],
                              stats    = ranked_genes_lst,
                              eps      = 0.0,
                              minSize  = 10,
                              maxSize  = 1500)
}

outF <- paste(outFolder, "/sleuth.results.gsea.pa.txt", sep="")
fwrite(fgseaRes$pa, file = outF, sep="\t", sep2=c("", ",", ""))

