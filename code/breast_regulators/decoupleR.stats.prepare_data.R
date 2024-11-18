
## metadata of the full seurat object

coi <- read_tsv("../breast_cell-cycle/coi.short.txt", col_names = FALSE) %>% pull(1)
coi_all <- read_tsv("../breast_cell-cycle/coi.long.txt", col_names = FALSE) %>% pull(1)

meta <- read_tsv("../breast_cell-cycle/so_metadata_post_integration.tsv.gz")

#fix colnames (strip 1s)
colnames(meta) <- gsub("1$", "", colnames(meta))

#filter for coi_all
meta <- meta %>% dplyr::filter(cell_type %in% coi_all)

#re-order factor levels of cell_type
meta <- meta %>% mutate(cell_type = factor(cell_type, levels = coi_all))

#strip MCF7
meta <- meta %>% mutate(cell_type_2 = gsub("MCF7_|MCF7-", "", cell_type))
meta <- meta %>% mutate(cell_type_2 = factor(cell_type_2, levels = gsub("MCF7_|MCF7-", "", coi_all)))

## metadata of the samples

metadata_ordered <- read.xls("../data/metadata_filtered_ordered.xlsx") %>% 
  as_tibble() %>%
  mutate(sample = paste0(description, "_batch-", batch)) %>%
  mutate(class = ifelse(ESR1_status == "mut", "ESR1_mutant", class)) %>%
  mutate(class = factor(class, levels = c("treatment-naive", "ESR1_mutant", "first_line", "resistant", "second_line"))) %>%
  dplyr::rename(cell_type = description) %>%
  mutate(cell_type_2 = gsub("MCF7_|MCF7-", "", cell_type))

## join the metadata(s)

meta_unique <- unique(metadata_ordered[,c("cell_type", "class", "treatment", "ESR1_status")])
meta_tmp <- merge(meta %>% rownames_to_column(var = "rowid"), meta_unique, by = "cell_type")
meta_tmp <- meta_tmp %>% column_to_rownames(var = "rowid")
meta_tmp <- meta_tmp[rownames(meta),]
meta$class <- meta_tmp$class
meta$treatment <- meta_tmp$treatment
meta$ESR1_status <- meta_tmp$ESR1_status

## define groups of conditions

groups <- list()
groups[["all_first_line"]] <- c("RM", "WM-2d", "Lted", "Tam-2d", "TamR", "Fulv-2d", "FulvR", "MCF7-D538G_RM", "MCF7-D538G_WM-2d")
groups[["wm_second_line"]] <- c("RM", "WM-2d", "Lted", "Lted_Tam-7d", "Lted_Fulv-7d", "Lted_Palbo-7d", "Lted_CDK7i-7d")
groups[["tam_second_line"]] <- c("RM", "Tam-2d", "TamR", "TamR_WM-7d", "TamR_Fulv-7d", "TamR_Palbo-7d", "TamR_CDK7i-7d")

