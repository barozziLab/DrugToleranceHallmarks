
##

dataPrep <- function(id) {
  
  #load the specified dataset from cBioPortal and fix the OS/DFS survival
  
  d <- cBioDataPack(id, ask = FALSE)
  
  meta <- colData(d)
  
  if (length(grep("OS_MONTHS",names(meta))) > 0) {
    meta[meta$OS_MONTHS == "[Not Available]" & !is.na(meta$OS_MONTHS), "OS_MONTHS"] <- NA
    meta$OS_MONTHS <- as.numeric(meta$OS_MONTHS)
  }
  
  if (length(grep("DFS_MONTHS",names(meta))) > 0) {
    meta[meta$DFS_MONTHS == "[Not Available]" & !is.na(meta$DFS_MONTHS), "DFS_MONTHS"] <- NA
    meta$DFS_MONTHS <- as.numeric(meta$DFS_MONTHS)
  }
  
  colData(d) <- meta
  
  return(d)
  
}


##

signatureScoreCalc <- function(se, goi) {
  
  #calculate signature score as sum of the expression of the genes of interest (goi; entrez)
  #se must be a summarizedExperiment object, defined on entrez IDs
  
  #extract expression matrix
  em <- assay(se)
  
  #filter based on goi
  em_sig <- em[rownames(em) %in% goi,]
  
  #calculate the score
  sig_score <- colSums(em_sig, na.rm = TRUE)
  
  return(list(sig_score = sig_score, sig_expr = em_sig))
  
}


##

signatureScoreToMetadata <- function(d, sig_score) {
  
  #add signature score (sig_score; and higher/lower than median) to cBioPortal dataset's (d) metadata (colData)
  #return an updated cBioPortal dataset
  
  meta <- colData(d)
  
  scores <- as.data.frame(sig_score)[rownames(meta),1]
  meta$signature_score <- scores
  
  score_median <- median(scores, na.rm = TRUE)
  meta$signature <- NA
  meta$signature[meta$signature_score >= score_median] <- "High"
  meta$signature[meta$signature_score < score_median] <- "Low"
  
  colData(d) <- meta
  
  return(d)
  
}


##

featureToMetadata <- function(d, featureValues, featureName) {

  meta <- colData(d)
  values <- as.data.frame(featureValues)[rownames(meta),1]
  meta[,featureName] <- values
  colData(d) <- meta
  
  return(d)
  
}


##

survRoutine <- function(d, survType) {
  
  #compute KM plot and stat for OS or DFS (survType)
  
  if (survType == "OS") {
    
    fit <- survfit(
      Surv(OS_MONTHS, as.numeric(substr(OS_STATUS, 1, 1))) ~ signature,
      data = colData(d)
    )
    test <- survdiff(
      Surv(OS_MONTHS, as.numeric(substr(OS_STATUS, 1, 1))) ~ signature,
      data = colData(d)
    )
    
  } else if (survType == "DFS") {
    
    fit <- survfit(
      Surv(DFS_MONTHS, as.numeric(substr(DFS_STATUS, 1, 1))) ~ signature,
      data = colData(d)
    )
    test <- survdiff(
      Surv(DFS_MONTHS, as.numeric(substr(DFS_STATUS, 1, 1))) ~ signature,
      data = colData(d)
    )
    
  } else if (survType == "RFS") {
    
    fit <- survfit(
      Surv(RFS_MONTHS, as.numeric(substr(RFS_STATUS, 1, 1))) ~ signature,
      data = colData(d)
    )
    test <- survdiff(
      Surv(RFS_MONTHS, as.numeric(substr(RFS_STATUS, 1, 1))) ~ signature,
      data = colData(d)
    )
    
  }
  
  return(list(data = colData(d), fit = fit, test = test))
  
}

