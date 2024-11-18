
library("edgeR")
library("gplots")

data <- read.table("GSE41245_processed_data.txt", sep="\t", header=T)
exp <- data[,3:ncol(data)]
gene_metadata <- data[,1:2]
rownames(exp) <- gene_metadata[,1]
sample_metadata <- read.table("GSE41245_metadata.txt", sep="\t", header=F)

#apply(exp, 2, sum)
#normalized to 1e6

treat_cols <- grep("EPCAM", colnames(exp))
ctrl_cols <- grep("IgG", colnames(exp))

#subtract each IgG from EPCAM
#good to calculate signature scores
#
exp_subtracted <- exp[,treat_cols] - exp[,ctrl_cols]
write.table(exp_subtracted, "GSE41245_processed_data.subtracted.txt", sep="\t", col.names=T, row.names=T, quote=F)
exp_logFC <- log2((exp[,treat_cols]+1) / (exp[,ctrl_cols]+1))
write.table(exp_logFC, "GSE41245_processed_data.log2FC.txt", sep="\t", col.names=T, row.names=T, quote=F)

#degs
#
degsPaired <- function(exp_data, design, logFC_t, p_t) {
    y <- DGEList(round(exp_data))
    y <- estimateDisp(y, design)
    y <- calcNormFactors(y, method="TMM")
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit)
    output <- topTags(qlf, 50000)$table
    degs_up_F <- output$PValue <= p_t & output$logFC >= logFC_t
    degs_up <- sort(rownames(output)[degs_up_F])
    degs_down_F <- output$PValue <= p_t & output$logFC <= -logFC_t
    degs_down <- sort(rownames(output)[degs_down_F])
    return(list(edger=output, degs_up=degs_up, degs_down=degs_down))
}
#
logFC_t <- log2(1.5)
p_t <- 0.05
#
#BCa, EPCAM vs IgG
#
f <- grep("cancer", sample_metadata[,4])
exp_data <- exp[,f]
exp_metadata <- sample_metadata[f,]
#
Subject <- factor(paste(as.numeric(exp_metadata[,3]), as.numeric(exp_metadata[,4]), sep="_"))
Treat <- unlist(strsplit(as.character(exp_metadata[,1]), split="\\."))
Treat <- factor(Treat[grep("row", Treat, invert=T)])
design <- model.matrix(~Subject+Treat)
#
out_cancer <- degsPaired(exp_data, design, logFC_t, p_t)
#up/down inverted
write.table(out_cancer$degs_up, file="DEGs.BCa.EPCAM_vs_IgG.down.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(out_cancer$degs_down, file="DEGs.BCa.EPCAM_vs_IgG.up.txt", sep="\t", quote=F, col.names=F, row.names=F)
#
#Healthy, EPCAM vs IgG
#
f <- grep("healthy", sample_metadata[,4])
exp_data <- exp[,f]
exp_metadata <- sample_metadata[f,]
#
Subject <- factor(paste(as.numeric(exp_metadata[,3]), as.numeric(exp_metadata[,4]), sep="_"))
Treat <- unlist(strsplit(as.character(exp_metadata[,1]), split="\\."))
Treat <- factor(Treat[grep("row", Treat, invert=T)])
design <- model.matrix(~Subject+Treat)
#
out_healthy <- degsPaired(exp_data, design, logFC_t, p_t)
write.table(out_healthy$degs_up, file="DEGs.Healthy.EPCAM_vs_IgG.down.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(out_healthy$degs_down, file="DEGs.Healthy.EPCAM_vs_IgG.up.txt", sep="\t", quote=F, col.names=F, row.names=F)
