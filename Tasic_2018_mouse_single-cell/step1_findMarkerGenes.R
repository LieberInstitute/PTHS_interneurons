library(BiocParallel)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Seurat)
library(future)
library(tidyverse)
library(writexl)

register(MulticoreParam(workers = parallel::detectCores()/2))
plan("multiprocess", workers = parallel::detectCores()/2)

load('rdas/Tasic2018_allen_mouse_cells_sce_sctranformed.rda')

########################################################################
##  marker genes by up-regulation in all pairwise cell type comparisons
markers.tTestUpAll = findMarkers(allen_mouse, test.type = 't', direction="up", 
                                 pval.type="all", groups = allen_mouse$subclass_label,
                                 row.data = rowData(allen_mouse))

###########################################################################
##  marker genes by up-regulation in 80% of pairwise cell type comparisons
markers.tTestUpSome = findMarkers(allen_mouse,  test.type = 't', direction="up", 
                                  pval.type="some", groups = allen_mouse$subclass_label, 
                                  lfc = 1, min.prop = .8, row.data = rowData(allen_mouse))

markerPvalbAll = markers.tTestUpAll[['Pvalb']] 
sapply(c(0.1, 0.05, 0.01), function(FDR) {
  rownames(markerPvalbAll)[markerPvalbAll$FDR <= FDR]
})
sapply(c(0.1, 0.05, 0.01), function(FDR) sum(markerPvalbAll$FDR <= FDR))
sapply(c(0.1, 0.05, 0.01), function(FDR) {
  grep('Kcnip2',rownames(markerPvalbAll)[markerPvalbAll$FDR <= FDR], value = T)
})


markerPvalbSome = markers.tTestUpSome[['Pvalb']] 
sapply(c(0.1, 0.05, 0.01), function(FDR) sum(markerPvalbSome$FDR <= FDR))
sapply(c(0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001), function(FDR) {
  grep('Kcnip2',rownames(markerPvalbSome)[markerPvalbSome$FDR <= FDR], value = T)
})

save(markers.tTestUpAll, markers.tTestUpSome, file = 'rdas/Tasic2018_allen_mouse_markerGenes.rda')


########################
# export files to excel
# filter markers by FDR < 0.05
alpha = 0.05
markers.All = lapply(markers.tTestUpAll, function(x) as.data.frame(x) %>% rownames_to_column(var = "Gene") %>% filter(FDR < alpha))

xlsxFileName = "tables/AllenMouse_markerGenes_UpInAllPairwise.xlsx"
write_xlsx(markers.All, path = xlsxFileName)
