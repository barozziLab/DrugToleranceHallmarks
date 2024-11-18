## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: decoupleR.stats.main.R
##
## Description: 
#~
## Wrapper for the regulator analysis
##
## Authors: 
#~
## Stephan Gruener & Iros Barozzi
##
## License: 
#~
## GNU GPL v3
## Copyright 2024 
## Copyright Iros Barozzi Stephan Gruener
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Notes:
#~
##
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(tidyverse)
library(gdata)
library(Seurat)
library(decoupleR)
library(msigdbr)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)

source("../lib/inc.func.decoupleR.R")


#~~~~~~~~~~~~~~~~~
#Prepare Gene Sets
##################

source("decoupleR.stats.prepare_gs.R")


#~~~~~~~~~~~~~~~~
#Data Preparation
#################

source("decoupleR.stats.prepare_data.R")

#metadata (comparisons)
metadata_comparison <- read.xls("../data/metadata_comparisons.xlsx") %>%
  as_tibble() %>%
  mutate(treated = gsub("MCF7_|MCF7-", "", treated)) %>%
  mutate(untreated = gsub("MCF7_|MCF7-", "", untreated))

# Thresholds
filt_tf_expr <- 0.05
filt_tf_act  <- 0.05


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Data Preparation - Neo-adjuvant
################################

# TF activities predictions

d_neo <- read_tsv("data_neo/neo_novel.addition.decoupleR.tfs.txt")
d_neo <- d_neo %>% select(-statistic) %>% dplyr::rename(gene = source)

# Feature importance from the models

d_neo_ml_feat_lasso <- read_tsv("data_neo/neo_novel.addition.ML.lasso.features.importance.txt")
d_neo_ml_feat_rf <- read_tsv("data_neo/neo_novel.addition.ML.rf.features.importance.txt")

d_neo <- d_neo %>% 
  left_join(d_neo_ml_feat_lasso %>% select(term, estimate) %>% dplyr::rename(gene = term, lasso_est = estimate), by = "gene") %>%
  left_join(d_neo_ml_feat_rf %>% dplyr::rename(gene = Variable, rf_imp = Importance), by = "gene")


#~~~~~~~~~~~~~~~~~~~~~~
#decoupleR Calculations
#######################

## Signaling

sig_all <- readRDS("decoupleR.activities.sig.rds")
sig_all <- sig_all %>%
  rowwise() %>%
  mutate(cell_type = strsplit(condition, split = "_batch")[[1]][1]) %>%
  ungroup() %>%
  mutate(cell_type = gsub("MCF7_|MCF7-", "", cell_type)) %>%
  mutate(cell_type = factor(cell_type, levels = gsub("MCF7_|MCF7-", "", coi_all))) %>%
  dplyr::filter(!is.na(cell_type))

## TFs

#run as part of decoupleR.calculations.tfs.vm.step2.R
tf_stats <- readRDS("decoupleR.activities.tfs.merge.tf_stats_merge.rds")

#select and relevel cell types
tf_stats <- tf_stats %>% 
  dplyr::filter(cell_type %in% coi_all) %>% 
  mutate(cell_type = gsub("MCF7_|MCF7-", "", cell_type)) %>%
  mutate(cell_type = factor(cell_type, levels = gsub("MCF7_|MCF7-", "", coi_all)))


#~~~~~~~~~~~~~
#Expressed TFs
##############

## Pre-processing

#run as part of decoupleR.calculations.tfs.vm.step2.R
tf_expr_frac <- readRDS("decoupleR.activities.tfs.merge.tf_expr_frac.rds")

#select and relevel cell types
tf_expr_frac <- tf_expr_frac %>% 
  dplyr::filter(cell_type %in% coi_all) %>% 
  mutate(cell_type = gsub("MCF7_|MCF7-", "", cell_type)) %>%
  mutate(cell_type = factor(cell_type, levels = gsub("MCF7_|MCF7-", "", coi_all)))

## General filter for expressed TFs in the dataset

# Filter for TFs expressed
tf_w_expr <- tf_expr_frac %>% 
  dplyr::filter(fraction >= filt_tf_expr) %>% 
  pull(TF_gene)

# Filter for TFs being significantly active
tf_w_act <- tf_stats %>% 
  dplyr::filter(n_cell_frac >= filt_tf_act) %>% 
  pull(source)

tf_w <- intersect(tf_w_expr, tf_w_act) %>% sort() %>% unique()


#~~~~~~~~~~~~~~~~~~~~~~
#TFs - Data Preparation
#######################

#stats: mean TF activity per (TF, condition)
tf_mat_mean <- tf_stats %>% 
  select(source, cell_type, mean) %>% 
  pivot_wider(names_from = "cell_type", values_from = "mean")

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
d_sub <- d_sub[, colnames(d_sub) != "RM"]
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
  mutate(name = factor(name, levels = gsub("MCF7_|MCF7-", "", coi_all)))

## Delta vs RM

d <- tf_mat_mean %>% column_to_rownames("source")
d_sub_rm <- d - d[,"RM"]
d_sub_rm <- d_sub_rm[, colnames(d_sub_rm) != "RM"]

# To long format
d_sub_rm_long <- d_sub_rm %>% 
  rownames_to_column(var = "gene") %>% 
  as_tibble() %>% 
  pivot_longer(cols = -c(gene))

# Re-order names
d_sub_rm_long <- d_sub_rm_long %>% 
  mutate(name = factor(name, levels = gsub("MCF7_|MCF7-", "", coi_all)))


#~~~~~~~~~~~~~~~~~~~
#TFs - Visualization
####################

# Plots based on WM 2d differences

source("decoupleR.stats.run_subtracted-paired.R")

# Plots based on all differences

source("decoupleR.stats.run_subtracted-rm.R")

# Seurat plots of selected markers/regulators

source("decoupleR.stats.run_markers.R")

