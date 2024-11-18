## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: run_downstream.R
##
## Description: 
#~
## Downstream analyses and visualizations for celltype predictions in breast data using TOSICA
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

source("lib.R")

source("celltype_predictions_viz.R")

if (!dir.exists("cell_type_vs_state_switch")) {dir.create("cell_type_vs_state_switch")}
source("cell_type_vs_cell_state.R")