#!/bin/bash

#Collect predicted statistics (VM/HPC)
R --vanilla < decoupleR.calculations.tfs.vm.step0.R
R --vanilla < decoupleR.calculations.tfs.vm.step1.R
R --vanilla < decoupleR.calculations.tfs.vm.step2.R

#Run downstream analyses (local)
R --vanilla < decoupleR.stats.main.R
