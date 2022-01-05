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

# get marker genes
alpha = 1
markers.All = lapply(markers.tTestUpAll, function(x) as.data.frame(x) %>% rownames_to_column(var = "Gene") %>% filter(FDR < alpha))
indKeep.AllMarkers = which(sapply(markers.All, nrow) > 0)
markerAll = unique(unlist(lapply(markers.All[indKeep.AllMarkers],'[[', 'Gene')))

indList = split(seq(ncol(allen_mouse)), allen_mouse$subclass_label)
counts = assays(allen_mouse)$logcounts
counts = counts[markerAll,]

# take average log10 gene levels by celltype
MF = do.call('cbind',lapply(indList, function(i){
  rowMeans(counts[,i])
}))

###################################
# load in PTHS mouse mega DESeq object #
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mega_tcf4_ages_DESeq2_svaAdj.rda')
load('/users/bphan/tcf4/PTHS_mouse/tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda')

# sample gene expression from PTHS mouse RNA-seq
geneVstAdult = varianceStabilizingTransformation(geneDds$Adult, blind=FALSE) # in case variances different by design
geneVstP1 = varianceStabilizingTransformation(geneDds$p1, blind=FALSE) # in case variances different by design
rownames(geneVstAdult) = outGeneList[[1]][rownames(geneVstAdult),'Symbol']
rownames(geneVstP1) = outGeneList[[1]][rownames(geneVstP1),'Symbol']
shared_genes = intersect(rownames(geneVstAdult), rownames(geneVstP1))
sample <-  cbind(assays(geneVstP1[shared_genes,])[[1]], 
                 assays(geneVstAdult[shared_genes,])[[1]]) # VST in log2
sample = sample / log2(10) # convert to log10

###############
# shared genes
shared_genes2 = intersect(rownames(sample), rownames(MF))
im_ref =as.matrix( MF[shared_genes2, ])
sample_out = sample[shared_genes2, ] 

write.csv(im_ref,"tables/Allen_Mouse_Brain_Reference_logCounts_p1-Adult.csv")
write.csv(sample_out,"tables/PTHSmouse_samples_logCounts_p1-Adult.csv")
