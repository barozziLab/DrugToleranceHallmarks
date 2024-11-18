
library("gdata")
library("psych")
library("MAST")
library("ROCR")

sample_metadata <- read.table("GSE51827_metadata.txt", sep="\t", header=F, comment.char="%")

data <- read.xls("GSE51827_readCounts.xls", sheet = 1, header = TRUE)
gene_metadata <- read.xls("GSE51827_platform.xls", sheet = 1, header = TRUE)

#sum by gene symbol
data_tmp <- merge(data, gene_metadata[,c("ID", "symbol")], by.x=1, by.y=1)
data_tmp <- data_tmp[data_tmp$symbol != "",]
exp <- data_tmp[,2:(ncol(data_tmp)-1)]
symbol <- as.character(data_tmp[,ncol(data_tmp)])
data_tmp_sum <- c()
for (i in 1:ncol(exp)) {
    tmp <- describeBy(exp[,i], group=symbol, mat=T)
    data_tmp_sum <- cbind(data_tmp_sum, tmp$n * tmp$mean)
}
rownames(data_tmp_sum) <- tmp$group1
colnames(data_tmp_sum) <- sample_metadata[,1]

#filter zeros
data_final <- data_tmp_sum[apply(data_tmp_sum, 1, sum) > 0, ]

#some cells are very lowly sequenced
#barplot(apply(data_final, 2, sum))

#avoid normalization
#normalize as needed [e.g. after computing pathway scores]
write.table(data_final, file="GSE51827_data.txt", sep="\t", quote=F, col.names=T, row.names=T)

#DEGs using MAST
#
degsUp <- function(expMat, cellsPos, cellsNeg, lbl_pos, lbl_neg, outPrefixMAST, foi) {
    #
    group <- colnames(expMat)
    group[cellsPos] <- lbl_pos
    group[cellsNeg] <- lbl_neg
    group <- group[cellsPos | cellsNeg]
    #
    expMat_subset <- expMat[,cellsPos | cellsNeg]
    #
    #filter out genes detected in less than 10 cells
    expMat_geneSum <- apply(expMat_subset > 0, 1, sum)
    expMat_MAST <- expMat_subset[expMat_geneSum > 10,]
    #
    #normalize by seq depth
    expMat_MAST <- t(t(expMat_MAST) / colSums(expMat_MAST))
    #
    raw <- FromMatrix(exprsArray=expMat_MAST, cData=data.frame(wellKey=colnames(expMat_MAST), group=group), fData=data.frame(primerid=rownames(expMat_MAST)))
    #
    two.sample <- LRT(raw, "group", referent=lbl_neg)
    p.bh <- p.adjust(two.sample[,"p.value"], method="BH")
    #
    #add fraction of cells expressing the gene (both classes)
    count <- apply(expMat_MAST[,group == lbl_pos] > 0, 1, sum)
    frac_pos <- count / length(group[group == lbl_pos])
    frac_pos <- frac_pos[order(names(frac_pos))]
    count <- apply(expMat_MAST[,group == lbl_neg] > 0, 1, sum)
    frac_rest <- count / length(group[group == lbl_neg])
    frac_rest <- frac_rest[order(names(frac_rest))]
    #
    #AUC
    pred_group <- rep(0, length(group))
    pred_group[group == lbl_pos] <- 1
    pred_auc <- apply(expMat_MAST, 1, auc, pred_group)
    pred_auc <- pred_auc[order(names(pred_auc))]
    #
    two.sample <- data.frame(two.sample, FDR=p.bh, expressingFracPos=round(frac_pos, 3), expressingFracNeg=round(frac_rest, 3), AUC=pred_auc)
    two.sample <- two.sample[order(two.sample$FDR),]
    #
    keep <- two.sample[,"direction"] == 1 & two.sample[,"p.value"] <= 0.05
    outPath <- paste(outPrefixMAST, ".txt", sep="")
    #
    out <- two.sample[keep, c(foi, "AUC", "expressingFracPos", "expressingFracNeg")]
    keep <- out[,"AUC"] >= 0.6
    out <- out[keep,]
    write.table(out, file=outPath, sep="\t", quote=F, row.names=F)
    return(out)
}
#
coi <- c("primerid", "lrstat", "p.value", "FDR")
#
classPos <- "CL"
classNeg <- "SC"
cellsPos <- sample_metadata[,4] == classPos
cellsNeg <- sample_metadata[,4] == classNeg
#
outPrefixMAST <- paste("DEGs.", classPos, "_vs_", classNeg, ".up", sep="")
out <- degsUp(data_final, cellsPos, cellsNeg, classPos, classNeg, outPrefixMAST, coi)
#
outPrefixMAST <- paste("DEGs.", classNeg, "_vs_", classPos, ".up", sep="")
out <- degsUp(data_final, cellsNeg, cellsPos, classNeg, classPos, outPrefixMAST, coi)
