
library("psych")
library("stringr")
library("limma")

#cleanup metadata
#
metadata <- read.csv("GSE20181_series_matrix.metadata.samples.txt", sep="\t", header=F)
patientID <- gsub("A|B|C", "", metadata[,2])
treat <- as.character(metadata[,5])
treat[treat == "Letrozole, 2.5mg/day,oral, 10-14 days"] <- "treatment_2w"
treat[treat == " Letrozole, 2.5mg/day,oral, 90 days"] <- "treatment_3mo"
response <- as.character(metadata[,6])
response[response == " responder"] <- "responder"
response[response == " nonresponder"] <- "nonresponder"
response[response == " not assessable"] <- "not_assessable"
metadata <- data.frame(ID=metadata[,1], patientID=patientID, treatment=treat, response=response)
write.table(metadata, file="GSE20181_series_matrix.metadata.samples.clean.txt", sep="\t", quote=F, col.names=T, row.names=F)

#annotate data table with gene symbol
#
data <- read.csv("GSE20181_series_matrix.data.txt", sep="\t", header=T)
plat_meta <- read.csv("../GSE5462_AI_early_response/GPL96-57554.txt", sep="\t", comment.char="#")
#
data_anno <- merge(data, plat_meta[,c(1,11)], by.x="ID_REF", by.y="ID")
#
#explode genes
out <- c()
for (i in 1:nrow(data_anno)) {
    gene_split <- strsplit(as.character(data_anno[i,"Gene.Symbol"]), " /// ")[[1]]
    if (length(gene_split) > 0) {
        for (j in 1:length(gene_split)) {
            out <- rbind(out, cbind(data_anno[i,], gene_split[j]))
        }
    }
}
#
#log2
exp <- out[,2:177]
exp <- log2(exp+1)
#
#average genes
tmp <- describeBy(exp, group=as.character(out[,ncol(out)]))
exp_avg <- c()
for (name in names(tmp)) exp_avg <- rbind(exp_avg, tmp[name][[1]][,"mean"])
rownames(exp_avg) <- names(tmp)
colnames(exp_avg) <- rownames(tmp[name][[1]])
#
write.table(round(exp_avg,2), file="GSE20181_series_matrix.data.clean.txt", sep="\t", quote=F, col.names=T, row.names=T)
#
#DEGs
#
#threshold on uncorrect p-value
#
runLimma <- function(expr, design, contrast, logFC_t, p_t) {
    fit <- lmFit(expr, design)
    fit <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit)
    output <- topTable(fit2, coef = 1, n=50000)
    degs_up_F <- output$P.Value <= p_t & output$logFC >= logFC_t
    degs_up <- sort(rownames(output)[degs_up_F])
    degs_down_F <- output$P.Value <= p_t & output$logFC <= -logFC_t
    degs_down <- sort(rownames(output)[degs_down_F])
    return(list(limma=output, degs_up=degs_up, degs_down=degs_down))
}
#
runLimmaPaired <- function(expr, classes, class_, sample_, logFC_t, p_t) {
    filter <- class_ %in% classes
    expr <- expr[,filter]
    sample_ <- as.factor(as.character(sample_[filter]))
    class_ <- as.factor(as.character(class_[filter]))
    design <- model.matrix(~ sample_ + class_)
    fit <- lmFit(expr, design)
    fit2 <- eBayes(fit)
    coef_id <- paste("class_", classes[1], sep="")
    output <- topTable(fit2, coef = coef_id, n=50000)
    degs_up_F <- output$P.Value <= p_t & output$logFC >= logFC_t
    degs_up <- sort(rownames(output)[degs_up_F])
    degs_down_F <- output$P.Value <= p_t & output$logFC <= -logFC_t
    degs_down <- sort(rownames(output)[degs_down_F])
    return(list(limma=output, degs_up=degs_up, degs_down=degs_down))
}
#
logFC_t <- log2(1.5)
p_t <- 0.05
#
class_ <- as.factor(paste(metadata[,3], metadata[,4], sep="_"))
sample_ <- as.factor(metadata[,2])
#
#UNPAIRED
design <- model.matrix(~0 + class_)
#
#responder vs non responder [pre-treatment]
contrast <- makeContrasts(class_pretreatment_responder - class_pretreatment_nonresponder, levels = design)
limmaOut_R_vs_NR <- runLimma(exp_avg, design, contrast, logFC_t, p_t)
write.table(limmaOut_R_vs_NR$limma, file="DEGs.Resp_Vs_NonResp.pretreat.limma.txt", sep="\t", quote=F, col.names=T, row.names=T)
write.table(limmaOut_R_vs_NR$degs_up, file="DEGs.Resp_Vs_NonResp.pretreat.up.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(limmaOut_R_vs_NR$degs_down, file="DEGs.Resp_Vs_NonResp.pretreat.down.txt", sep="\t", quote=F, col.names=F, row.names=F)
#
#responder vs non responder [2w]
contrast <- makeContrasts(class_treatment_2w_responder - class_treatment_2w_nonresponder, levels = design)
limmaOut_R_vs_NR <- runLimma(exp_avg, design, contrast, logFC_t, p_t)
write.table(limmaOut_R_vs_NR$limma, file="DEGs.Resp_Vs_NonResp.treat2w.limma.txt", sep="\t", quote=F, col.names=T, row.names=T)
write.table(limmaOut_R_vs_NR$degs_up, file="DEGs.Resp_Vs_NonResp.treat2w.up.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(limmaOut_R_vs_NR$degs_down, file="DEGs.Resp_Vs_NonResp.treat2w.down.txt", sep="\t", quote=F, col.names=F, row.names=F)
#
#responder vs non responder [3mo]
contrast <- makeContrasts(class_treatment_3mo_responder - class_treatment_3mo_nonresponder, levels = design)
limmaOut_R_vs_NR <- runLimma(exp_avg, design, contrast, logFC_t, p_t)
write.table(limmaOut_R_vs_NR$limma, file="DEGs.Resp_Vs_NonResp.treat3mo.limma.txt", sep="\t", quote=F, col.names=T, row.names=T)
write.table(limmaOut_R_vs_NR$degs_up, file="DEGs.Resp_Vs_NonResp.treat3mo.up.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(limmaOut_R_vs_NR$degs_down, file="DEGs.Resp_Vs_NonResp.treat3mo.down.txt", sep="\t", quote=F, col.names=F, row.names=F)
#
#PAIRED
#
#responder [2w]
classes <- c("treatment_2w_responder", "pretreatment_responder")
limmaOut_T_vs_P <- runLimmaPaired(exp_avg, classes, class_, sample_, logFC_t, p_t)
write.table(limmaOut_T_vs_P$limma, file="DEGs.Resp.Treat_Vs_Pretreat.treat2w.limma.txt", sep="\t", quote=F, col.names=T, row.names=T)
write.table(limmaOut_T_vs_P$degs_up, file="DEGs.Resp.Treat_Vs_Pretreat.treat2w.up.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(limmaOut_T_vs_P$degs_down, file="DEGs.Resp.Treat_Vs_Pretreat.treat2w.down.txt", sep="\t", quote=F, col.names=F, row.names=F)
#responder [3mo]
classes <- c("treatment_3mo_responder", "pretreatment_responder")
limmaOut_T_vs_P <- runLimmaPaired(exp_avg, classes, class_, sample_, logFC_t, p_t)
write.table(limmaOut_T_vs_P$limma, file="DEGs.Resp.Treat_Vs_Pretreat.treat3mo.limma.txt", sep="\t", quote=F, col.names=T, row.names=T)
write.table(limmaOut_T_vs_P$degs_up, file="DEGs.Resp.Treat_Vs_Pretreat.treat3mo.up.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(limmaOut_T_vs_P$degs_down, file="DEGs.Resp.Treat_Vs_Pretreat.treat3mo.down.txt", sep="\t", quote=F, col.names=F, row.names=F)
#
#non responder [2w]
classes <- c("treatment_2w_nonresponder", "pretreatment_nonresponder")
limmaOut_T_vs_P <- runLimmaPaired(exp_avg, classes, class_, sample_, logFC_t, p_t)
write.table(limmaOut_T_vs_P$limma, file="DEGs.NonResp.Treat_Vs_Pretreat.treat2w.limma.txt", sep="\t", quote=F, col.names=T, row.names=T)
write.table(limmaOut_T_vs_P$degs_up, file="DEGs.NonResp.Treat_Vs_Pretreat.treat2w.up.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(limmaOut_T_vs_P$degs_down, file="DEGs.NonResp.Treat_Vs_Pretreat.treat2w.down.txt", sep="\t", quote=F, col.names=F, row.names=F)
#non responder [3mo]
classes <- c("treatment_3mo_nonresponder", "pretreatment_nonresponder")
limmaOut_T_vs_P <- runLimmaPaired(exp_avg, classes, class_, sample_, logFC_t, p_t)
write.table(limmaOut_T_vs_P$limma, file="DEGs.NonResp.Treat_Vs_Pretreat.treat3mo.limma.txt", sep="\t", quote=F, col.names=T, row.names=T)
write.table(limmaOut_T_vs_P$degs_up, file="DEGs.NonResp.Treat_Vs_Pretreat.treat3mo.up.txt", sep="\t", quote=F, col.names=F, row.names=F)
write.table(limmaOut_T_vs_P$degs_down, file="DEGs.NonResp.Treat_Vs_Pretreat.treat3mo.down.txt", sep="\t", quote=F, col.names=F, row.names=F)
