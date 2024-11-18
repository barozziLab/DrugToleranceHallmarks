#!/bin/bash

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## Script name: run_tosica.sh
##
## Description: 
#~
## Wrapper for cell type predictions
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
## install TOSICA environment before running: conda env create --name TOSICA --file=TOSICA.yml
## run on HPC, preferrably on a GPU partition
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

conda activate TOSICA

python train_predict_full_epithelial.py
python train_predict_broad.py