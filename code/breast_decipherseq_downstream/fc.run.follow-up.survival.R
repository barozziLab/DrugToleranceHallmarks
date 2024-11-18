
options(connectionObserver = NULL)

library("org.Hs.eg.db")
library("cBioPortalData")
library("survival")
library("survminer")
library("cowplot")

source("../lib/inc.func.survival.R")

out_folder_surv <- "fc_results_survival/"
cmd <- paste0("mkdir ", out_folder_surv)
system(cmd)

#plot titles
plt_titles <- list()
plt_titles[["OS"]] <- "OS"
plt_titles[["RFS"]] <- "RFS"
plt_titles[["DFS"]] <- "DFS"


#~~~~~~~~~
#Signature
##########

mrk_list_merge <- mrk_list

#PA
for (sig_id in c("PA_UP", "PA_DOWN")) {
  mrk_list_merge[[sig_id]] <- m_t2g[["PA"]] %>% 
    as_tibble() %>% 
    dplyr::filter(gs_name == sig_id) %>%
    pull(gene_symbol)
}

#IFN hallmark
for (sig_id in c("HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE")) {
  sig_id_short <- gsub("HALLMARK_", "", sig_id)
  mrk_list_merge[[sig_id_short]] <- m_t2g[["Hallmarks"]] %>% 
    as_tibble() %>% 
    dplyr::filter(gs_name == sig_id) %>%
    pull(gene_symbol)
}

#ISGF3
mrk_list_merge[["ISGF3"]] <- c("IRF9", "STAT1", "STAT2")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#cBioPortal Data Preparation
############################

#see: https://waldronlab.io/cBioPortalData/articles/cBioPortalData.html

#to clean cached data use: 
#unlink("~/.cache/cBioPortalData/")

cbio <- cBioPortal(hostname = "www.cbioportal.org", protocol = "https", api. = "/api/v2/api-docs")
studies <- getStudies(cbio)

#cancer types of potential interest --> cancerTypeId
#breast cancer                      --> brca

#studies %>% filter(cancerTypeId == "brca")

#only consider the largest
#i.e. c("brca_metabric", "brca_tcga")


#~~~~~~~~~~~~~~~~~~~~~~~~~
#cBioPortal: brca_metabric
##########################

study_id <- "brca_metabric"
rna_slot <- "mrna_illumina_microarray_zscores_ref_diploid_samples"

d_all <- dataPrep(study_id)

#filter on ER/HER2 status + hormone therapy
#keep only TNBC
w <- d_all$ER_STATUS == "Positive" & d_all$PR_STATUS == "Positive" & d_all$HER2_STATUS == "Negative" & d_all$HORMONE_THERAPY == "YES"
d_lum <- d_all[,w]

#endpoint max 120 months
thresh <- 120
meta <- colData(d_lum)
meta_mod <- colData(d_lum)
meta_mod$OS_MONTHS[meta$OS_MONTHS > thresh] <- thresh
meta_mod$OS_STATUS[meta$OS_MONTHS > thresh] <- "0:LIVING"
meta_mod$RFS_MONTHS[meta$RFS_MONTHS > thresh] <- thresh
meta_mod$RFS_STATUS[meta$RFS_STATUS > thresh] <- "0:LIVING"
colData(d_lum) <- meta_mod

#extract mRNA expression data
se <- d_lum@ExperimentList@listData[[rna_slot]]

#survival

survRes <- list()
plts    <- list()

summary_stats <- c()

for (sig_id in names(mrk_list_merge)) {
  
  survRes[[sig_id]] <- list()
  plts[[sig_id]] <- list()
  
  #calculate and store the signature scores
  sig_res <- signatureScoreCalc(se, mrk_list_merge[[sig_id]])
  d_lum <- signatureScoreToMetadata(d_lum, sig_res$sig_score)
  
  for (surv_type in c("OS", "RFS")) {
  
    survRes[[sig_id]][[surv_type]] <- survRoutine(d_lum, surv_type)
    
    plts[[sig_id]][[surv_type]] <- ggsurvplot(survRes[[sig_id]][[surv_type]]$fit, 
                                              data = survRes[[sig_id]][[surv_type]]$data, 
                                              risk.table = TRUE, 
                                              pval = TRUE, 
                                              pval.method= TRUE,
                                              pval.size = 3,
                                              pval.method.size = 3,
                                              ggtheme = theme_bw(), 
                                              title = paste0("Program ", sig_id, " - ", plt_titles[[surv_type]]))
    
    test_res <- survRes[[sig_id]][[surv_type]]$test
    pval <- pchisq(test_res$chisq, length(test_res$n)-1, lower.tail = FALSE)
    summary_stats <- rbind(summary_stats, c(sig_id, surv_type, test_res$chisq, pval))
  
  }
  
}

summary_stats_tib <- tibble(program = summary_stats[,1],
                            surv_type = summary_stats[,2],
                            chisq = as.numeric(summary_stats[,3]),
                            p = as.numeric(summary_stats[,4]))

summary_stats_tib$p.adjust <- p.adjust(summary_stats_tib$p, method = "BH")

#save summary
outF <- paste0(out_folder_surv, "/brca_metabric.programs.summary.txt", sep = "")
write_tsv(x = summary_stats_tib, file = outF)

#save plots of the genes with a significant difference (super-lenient)
programs_plt <- summary_stats_tib %>% 
  dplyr::filter(p <= 0.25)

for (surv_type_i in c("OS", "RFS")) {

  programs_plt_f <- programs_plt %>% 
    dplyr::filter(surv_type == surv_type_i)
  outF <- paste0(out_folder_surv, "/brca_metabric.programs.km.", surv_type_i, ".top.pdf", sep = "")
  pdf(outF, width = 6, height = 5)
  for (i in 1:nrow(programs_plt_f)) {
    poi <- programs_plt_f[i,] %>% pull(program)
    soi <- programs_plt_f[i,] %>% pull(surv_type)
    print(plts[[poi]][[soi]])
  }
  dev.off()

}


#~~~~~~~~~~~~~~~~~~~~~
#cBioPortal: brca_tcga
######################

study_id <- "brca_tcga"
rna_slot <- "mrna_seq_v2_rsem_zscores_ref_diploid_samples"

d_all <- dataPrep(study_id)

#filter on ER/HER2 status only
#no filter on PHARMACEUTICAL_TX_ADJUVANT, numbers too low
#keep only TNBC
w <- d_all$ER_STATUS_BY_IHC == "Positive" & d_all$IHC_HER2 == "Positive"
d_lum <- d_all[,w]

#extract mRNA expression data
se <- d_lum@ExperimentList@listData[[rna_slot]]

#endpoint max 120 months
thresh <- 120
meta <- colData(d_lum)
meta_mod <- colData(d_lum)
meta_mod$OS_MONTHS[meta$OS_MONTHS > thresh] <- thresh
meta_mod$OS_STATUS[meta$OS_MONTHS > thresh] <- "0:LIVING"
meta_mod$DFS_MONTHS[meta$DFS_MONTHS > thresh] <- thresh
meta_mod$DFS_STATUS[meta$DFS_STATUS > thresh] <- "0:LIVING"
colData(d_lum) <- meta_mod

#extract mRNA expression data
se <- d_lum@ExperimentList@listData[[rna_slot]]

#survival

survRes <- list()
plts    <- list()

summary_stats <- c()

for (sig_id in names(mrk_list_merge)) {
  
  survRes[[sig_id]] <- list()
  plts[[sig_id]] <- list()
  
  #calculate and store the signature scores
  sig_res <- signatureScoreCalc(se, mrk_list_merge[[sig_id]])
  d_lum <- signatureScoreToMetadata(d_lum, sig_res$sig_score)
  
  for (surv_type in c("OS", "DFS")) {
    
    survRes[[sig_id]][[surv_type]] <- survRoutine(d_lum, surv_type)
    
    plts[[sig_id]][[surv_type]] <- ggsurvplot(survRes[[sig_id]][[surv_type]]$fit, 
                                              data = survRes[[sig_id]][[surv_type]]$data, 
                                              risk.table = TRUE, 
                                              pval = TRUE, 
                                              pval.method= TRUE,
                                              pval.size = 3,
                                              pval.method.size = 3,
                                              ggtheme = theme_bw(), 
                                              title = paste0("Program ", sig_id, " - ", plt_titles[[surv_type]]))
    
    test_res <- survRes[[sig_id]][[surv_type]]$test
    pval <- pchisq(test_res$chisq, length(test_res$n)-1, lower.tail = FALSE)
    summary_stats <- rbind(summary_stats, c(sig_id, surv_type, test_res$chisq, pval))
    
  }
  
}

summary_stats_tib <- tibble(program = summary_stats[,1],
                            surv_type = summary_stats[,2],
                            chisq = as.numeric(summary_stats[,3]),
                            p = as.numeric(summary_stats[,4]))

summary_stats_tib$p.adjust <- p.adjust(summary_stats_tib$p, method = "BH")

#save summary
outF <- paste0(out_folder_surv, "/brca_tcga.programs.summary.txt", sep = "")
write_tsv(x = summary_stats_tib, file = outF)

#save plots of the genes with a significant difference (super-lenient)
programs_plt <- summary_stats_tib %>% 
  dplyr::filter(p <= 0.25)

for (surv_type_i in c("OS", "DFS")) {
  
  programs_plt_f <- programs_plt %>% 
    dplyr::filter(surv_type == surv_type_i)
  outF <- paste0(out_folder_surv, "/brca_tcga.programs.km.", surv_type_i, ".top.pdf", sep = "")
  pdf(outF, width = 6, height = 5)
  for (i in 1:nrow(programs_plt_f)) {
    poi <- programs_plt_f[i,] %>% pull(program)
    soi <- programs_plt_f[i,] %>% pull(surv_type)
    print(plts[[poi]][[soi]])
  }
  dev.off()
  
}

