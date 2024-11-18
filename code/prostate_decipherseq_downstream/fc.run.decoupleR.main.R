
library(Seurat)
library(decoupleR)
library(ggrepel)
library(gdata)

source("../lib/inc.func.decoupleR.R")


#~~~~~~~~~~~~~~~~
#Data Preparation
#################

so <- readRDS("so_sct_clusters_cc.rds")

# Change identity to cell type
Idents(so) <- so@meta.data$cell_type

#metadata
metadata_ordered <- read.xls("../data/metadata_filtered_ordered.xlsx", sheet = 2) %>% 
  as_tibble()

#metadata (comparisons)
metadata_comparison <- read.xls("../data/metadata_comparisons.xlsx") %>%
  as_tibble() %>%
  dplyr::filter(reduced == "Y") %>% 
  dplyr::select(treated, untreated, group)

#
cell_type_include <- c("LNCaP_RM", "LNCaP_WM-7d", "LNCaP95", "LNCaP95_WM-7d", "LNCaP95_Enz-7d")

# Thresholds
filt_tf_expr <- 0.05
filt_tf_act  <- 0.05


#~~~~~~~~~~~~~~~~~~~~~~
#decoupleR Calculations
#######################

acts_tf_path <- "decoupleR.activities.tfs.rds"
#acts_sig_path <- "decoupleR.activities.sig.rds"
so_decoupleR_path <- "so_sct_clusters_cc.decoupleR.rds"

if (file.exists(acts_tf_path)) {
  
  net_tf_acts <- readRDS(acts_tf_path)
  #net_sig_acts <- readRDS(acts_sig_path)  ## Frozen for memory issues
  so <- readRDS(so_decoupleR_path)
  
} else {
  
  source("fc.run.decoupleR.calculations.R")
  
}

# Re-sort cell type identifiers
so@meta.data$cell_type <- factor(so@meta.data$cell_type, levels = c("LNCaP_RM", "LNCaP_WM-7d", "LNCaP95", "LNCaP95_WM-7d", "LNCaP95_Enz-7d"))
Idents(so) <- so@meta.data$cell_type

mat <- as.matrix(so@assays$SCT@data)


#~~~~~~~~~~~~~
#Expressed TFs
##############

## Pre-processing

# Change default assay to tfsulm
DefaultAssay(object = so) <- "tfsulm"

# Scale the data
so <- ScaleData(so)
so@assays$tfsulm@data <- so@assays$tfsulm@scale.data

## Convert net_tf_acts to

# 1) fraction of cells per cell-type
# 2) average activity per cell-type

pvalue_thresh <- 0.05

meta_tmp <- so@meta.data %>%
  rownames_to_column(var = "cellID") %>% 
  select(cellID, cell_type) %>%
  as_tibble()

tfs_to_consider <- net_tf_acts %>% 
  filter(p_value <= pvalue_thresh) %>% 
  pull(source) %>% 
  unique()

tf_stats <- net_tf_acts %>% 
  mutate(significant = p_value <= pvalue_thresh) %>%
  dplyr::filter(source %in% tfs_to_consider) %>% 
  dplyr::rename(cellID = condition) %>%
  left_join(meta_tmp, by = "cellID") %>%
  group_by(source, cell_type) %>%
  summarise(mean = mean(score), n_cells_pos = sum(significant), n_cells_tot = n()) %>%
  mutate(n_cell_frac = n_cells_pos / n_cells_tot) %>% 
  ungroup()

out_f <- paste0(out_folder, "/decoupleR.tfs_stats.long.txt")
write_tsv(x = tf_stats, file = out_f)

#stats: mean TF activity per (TF, condition)
tf_mat_mean <- tf_stats %>% 
  select(source, cell_type, mean) %>% 
  pivot_wider(names_from = "cell_type", values_from = "mean")

#stats: fraction of cells with significant activity per (TF, condition)
tf_mat_frac <- tf_stats %>% 
  select(source, cell_type, n_cell_frac) %>% 
  pivot_wider(names_from = "cell_type", values_from = "n_cell_frac")

## Prepare tables with expression fractions

tf_expr_frac <- list()

tf_ids <- rownames(so@assays$tfsulm) %>% 
  unique() %>% 
  sort()

df_tf_gene <- t(mat[rownames(mat) %in% tf_ids,]) %>%
  as.data.frame() %>%
  mutate(cell_type = Idents(so)) %>%
  pivot_longer(cols = -c(cell_type), names_to = "TF_gene", values_to = "expression")

# Fraction of cells expressing the gene per cell type
tf_expr_frac[["cell-type"]] <- df_tf_gene %>%
  group_by(cell_type, TF_gene) %>%
  summarise(more0 = sum(expression > 0), tot = n()) %>%
  mutate(fraction = more0 / tot)

## General filter for expressed TFs in the datast

# Filter for TFs expressed
tf_w_expr <- tf_expr_frac$`cell-type` %>% 
  dplyr::filter(fraction >= filt_tf_expr) %>% 
  pull(TF_gene) %>%
  sort() %>% 
  unique()

# Filter for TFs being significantly active
tf_w_act <- tf_stats %>% 
  dplyr::filter(n_cell_frac >= filt_tf_act) %>% 
  pull(source) %>%
  sort() %>% 
  unique()

tf_w <- intersect(tf_w_expr, tf_w_act) %>% sort() %>% unique()


#~~~~~~~~~~~~~~~~~~~~~~
#TFs - Data Preparation
#######################

## Delta vs matched "base-line"

#see metadata_comparison

d <- tf_mat_mean %>% column_to_rownames("source")
d_sub <- d
for (i in 1:nrow(metadata_comparison)) {
  untreated_col <- metadata_comparison[i,] %>% pull(untreated)
  treated_col <- metadata_comparison[i,] %>% pull(treated)
  if (untreated_col %in% colnames(d) & treated_col %in% colnames(d)) {
    d_sub[,treated_col] <- d[,treated_col] - d[,untreated_col]
  }
}
d_sub <- d_sub[, colnames(d_sub) != "LNCaP_RM"]
cnames_sorted <- metadata_comparison$treated %>% unique()
cnames_sorted <- cnames_sorted[cnames_sorted %in% colnames(d)]
d_sub <- d_sub[,cnames_sorted]

# To long format
d_sub_long <- d_sub %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c(gene))

# Re-order names
d_sub_long <- d_sub_long %>% 
  mutate(name = factor(name, levels = cell_type_include[cell_type_include %in% d_sub_long$name]))

## Delta vs RM

d <- tf_mat_mean %>% column_to_rownames("source")
d_sub_rm <- d - d[,"LNCaP_RM"]
d_sub_rm <- d_sub_rm[, colnames(d_sub_rm) != "LNCaP_RM"]

# To long format
d_sub_rm_long <- d_sub_rm %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c(gene))

# Re-order names
d_sub_rm_long <- d_sub_rm_long %>% 
  mutate(name = factor(name, levels = cell_type_include[cell_type_include %in% d_sub_rm_long$name]))


#~~~~~~~~~~~~~~~~~~~
#TFs - Visualization
####################

# Plots based on WM 7d differences

source("fc.run.decoupleR.main.subtracted-paired.R")


