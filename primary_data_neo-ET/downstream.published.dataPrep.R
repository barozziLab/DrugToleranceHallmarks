
gs_neo <- list()

dirs <- c("DOWN", "UP")

#~~~~~~~~~~~~~~~~~~~~~~~~~~
#Prepare: HR+; Neo-adjuvant
###########################

inMainFolder <- "../literature_meta-analysis/Neo-Adjuvant/"

#GSE20181_AI_early_response

inSubFolder <- "GSE20181_AI_early_response"

#names_oi <- c("DEGs.Resp.Treat_Vs_Pretreat.treat2w", "DEGs.Resp.Treat_Vs_Pretreat.treat3mo")
names_oi <- c("DEGs.Resp.Treat_Vs_Pretreat.treat3mo")

#names_oi_short <- c("AI-2weeks_Responsive", "AI-3months_Responsive")
names_oi_short <- c("AI-3months_Responsive")

for (i in 1:length(names_oi)) {
	for (dir in dirs) {
		inF <- paste(inMainFolder, inSubFolder, "/", names_oi[i], ".", tolower(dir), ".txt", sep = "")
		sig <- read_tsv(inF, col_names = FALSE) %>% pull(1)
		name <- paste("HR_Neoadj_Miller_", names_oi_short[i], "_", dir, sep = "")
		gs_neo[[name]] <- sig
	}
}

#kastrati_et_al_tam_neoadjuvant

inSubFolder <- "kastrati_et_al_tam_neoadjuvant"

for (dir in dirs) {
	inF <- paste(inMainFolder, inSubFolder, "/kastrati_et_al.", tolower(dir), ".txt", sep = "")
	sig <- read_tsv(inF, col_names = FALSE) %>% pull(1)
	name <- paste("HR_Neoadj_Kastrati_TAM_", dir, sep = "")
	gs_neo[[name]] <- sig
}

#patani_et_al_neoadjuvant

inSubFolder <- "patani_et_al_neoadjuvant"

for (treat in c("AI", "FULV")) {
	for (dir in dirs) {
		inF <- paste(inMainFolder, inSubFolder, "/", treat, ".", tolower(dir), ".txt", sep = "")
		sig <- read_tsv(inF, col_names = FALSE) %>% pull(1)
		name <- paste("HR_Neoadj_Patani_", treat, "_", dir, sep = "")
		gs_neo[[name]] <- sig
	}
}

#turnbull_et_al_AI_neoadjuvant

inSubFolder <- "turnbull_et_al_AI_neoadjuvant"

d_turnbull <- paste(inMainFolder, inSubFolder, "/limma_results.rds", sep = "")

names_oi <- names(d_turnbull)
names_oi <- names_oi[2]
#names_oi_short <- c("AI-2weeks", "AI-3months")
names_oi_short <- c( "AI-3months")

for (i in 1:length(names_oi)) {
	name_oi <- names_oi[i]
	name <- paste("HR_Neoadj_Turnbull_", names_oi_short[i], "_DOWN", sep = "")
	gs_neo[[name]] <- d_turnbull[[name_oi]] %>% filter(logFC <= -1 & adj.P.Val <= 0.05 & Symbol != "") %>% arrange(Symbol) %>% pull(Symbol)
	name <- paste("HR_Neoadj_Turnbull_", names_oi_short[i], "_UP", sep = "")
	gs_neo[[name]] <- d_turnbull[[name_oi]] %>% filter(logFC >= 1 & adj.P.Val <= 0.05 & Symbol != "") %>% arrange(Symbol) %>% pull(Symbol)
}

###########################

#~~~~~~~~~~~~~~~
#Post-processing
################

#exclude empty sets
w <- which(t(sapply(gs_neo, length)) > 0)
gs_neo <- gs_neo[w]

#fix names
names(gs_neo) <- gsub("neoadj", "Neoadj", names(gs_neo))
names(gs_neo) <- gsub("-3months", "", names(gs_neo))

################
