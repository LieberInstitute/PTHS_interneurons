library(SingleCellExperiment)
library(scater)
library(scran)
library(BiocParallel)
library(Seurat)
library(future)
library(scrattch.io)

register(MulticoreParam(workers = parallel::detectCores()/2))
plan("multiprocess", workers = parallel::detectCores()/2)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 10*1024^3)

# read in matrix of gene expression (by exon mapping)
tome        <- "data/transcrip.tome"
exons <- read_tome_dgCMatrix(tome,"data/t_exon")    # (or data/exon)
rownames(exons)   <- read_tome_gene_names(tome)
colnames(exons) <- read_tome_sample_names(tome)  

# read in the sample and gene metadata
pd = read.csv('data/sample_annotations.csv')
rownames(pd) = pd$sample_name
pd = pd[colnames(exons),]

# read in the gene metadata
geneDat = read.csv('data/mouse_VISp_2018-06-14_genes-rows.csv')
rownames(geneDat) = geneDat$gene_symbol
geneDat = geneDat[rownames(exons),]

# create RSE experiment and cell annotations
allen_mouse = SingleCellExperiment(assays = list(counts = exons), colData = pd, rowData = geneDat)
allen_mouse = allen_mouse[,allen_mouse$class_label != 'Exclude'] # remove doublets/low-quality

save(allen_mouse, file = 'rdas/Tasic2018_allen_mouse_cells_sce.rda', compress = T)

# Scaling and Size-factor normalization w/ Seurat
allen_mouse_seurat <- as.Seurat(allen_mouse, data = 'counts')
allen_mouse_seurat = SCTransform(allen_mouse_seurat)
# allen_mouse_seurat <- NormalizeData(allen_mouse_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
# allen_mouse_seurat <- FindVariableFeatures(allen_mouse_seurat, selection.method = "vst", nfeatures = 2000)
# allen_mouse_seurat <- ScaleData(allen_mouse_seurat, features = rownames(allen_mouse_seurat))

# PCA reduce dimensions and perform UMAP embedding
allen_mouse_seurat <- RunPCA(allen_mouse_seurat, dims = 100, verbose = TRUE)
allen_mouse_seurat <- RunUMAP(allen_mouse_seurat, dims = 1:30, verbose = TRUE)
allen_mouse = as.SingleCellExperiment(allen_mouse_seurat)
save(allen_mouse_seurat, file = 'rdas/Tasic2018_allen_mouse_cells_sctranformed_seurat.rda', compress = T)
save(allen_mouse, file = 'rdas/Tasic2018_allen_mouse_cells_sctranformed_sce.rda', compress = T)


