library(DESeq2)
library(SingleCellExperiment)
library(tidyverse)
library(parallel)
library(future.apply)
library(future)
# reference gene expression profile is Allen scRNA-seq cells

ncores = parallel::detectCores()/2
plan("multiprocess", workers = parallel::detectCores()/2)
options(future.globals.maxSize = 10* 1024^3)

##################################################
# load in Tasic cell type specific marker genes #
load('rdas/Tasic2018_allen_mouse_cells_sce_sctranformed.rda')
load('rdas/Tasic2018_allen_mouse_markerGenes.rda')
