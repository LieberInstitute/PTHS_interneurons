library(DESeq2)
library(SingleCellExperiment)
library(tidyverse)
library(parallel)
library(future.apply)
library(future)
library(ggplot2)

# reference gene expression profile is Allen scRNA-seq cells

ncores = parallel::detectCores()/2
plan("multiprocess", workers = parallel::detectCores()/2)
options(future.globals.maxSize = 10* 1024^3)

##################################################
# load in Tasic cell type specific marker genes #
if(FALSE){
  load('rdas/Tasic2018_allen_mouse_cells_sce_sctranformed.rda')
  load('rdas/Tasic2018_allen_mouse_markerGenes.rda')
  
  # get Tcf4 gene expression
  cellData = as.data.frame(colData(allen_mouse))
  cellData$Tcf4 = assays(allen_mouse)$logcounts['Tcf4',]
  cellData = cellData[cellData$class_label != 'Exclude',]
  
  # change up label factors
  cellData$subclass_label = make.names(cellData$subclass_label)
  subclass_label = unique(cellData$subclass_label[order(cellData$subclass_order)])
  class_label = unique(cellData$class_label[order(cellData$class_order)])
  cellData$subclass_label = factor(cellData$subclass_label, levels = subclass_label)
  cellData$class_label = factor(cellData$class_label,  levels = class_label)
  
  # store for fast future loading
  save(cellData,file = 'rdas/Tasic2018_allen_mouse_Tcf4_expression.rda')
  
} else {
  load('rdas/Tasic2018_allen_mouse_Tcf4_expression.rda')
}



##################################################
# load in Tasic cell type specific marker genes #
subclass_color = cellData$subclass_color
names(subclass_color) = cellData$subclass_label
subclass_color = subclass_color[!duplicated(subclass_color)]

pdf('plots/Tasic2018_allen_mouse_Tcf4Expres_short.pdf', height = 3, width = 8.5)
ggplot(cellData, aes(x = subclass_label, y = 2^Tcf4, fill = subclass_label)) +
  geom_violin(position = 'dodge')+ #geom_point(position = position_jitterdodge(), pch = 21) +
  scale_fill_manual(values = subclass_color, guide = FALSE) + 
  facet_grid(~class_label, scales = 'free_x', space = 'free') + 
  theme_bw(base_size = 12) +   ylab("Normalized Tcf4 counts") + xlab('Allen Mouse Cell type') +
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1,size =8),
        axis.title.y = element_text(hjust = -.5))
dev.off()


pdf('plots/Tasic2018_allen_mouse_Tcf4Expres_wide.pdf', height = 5, width = 8)
ggplot(cellData, aes(x = subclass_label, y = 2^Tcf4, fill = subclass_label)) +
  geom_violin(position = 'dodge')+ #geom_point(position = position_jitterdodge(), pch = 21) +
  scale_fill_manual(values = subclass_color, guide = FALSE) + 
  facet_grid(~class_label, scales = 'free_x', space = 'free') + 
  theme_bw(base_size = 12) +   ylab("Normalized Tcf4 counts") + xlab('Allen Mouse Cell type') +
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1))
dev.off()




pdf('plots/Tasic2018_allen_mouse_Tcf4Expres_GABAergic.pdf', width = 3.5, height = 3)
ggplot(subset(cellData, class_label =='GABAergic'), 
       aes(x = subclass_label, y = 2^Tcf4, fill = subclass_label)) +
  geom_violin(position = 'dodge')+ #geom_point(position = position_jitterdodge(), pch = 21) +
  scale_fill_manual(values = subclass_color, guide = FALSE) + 
  #facet_grid(~class_label, scales = 'free_x', space = 'free') + 
  theme_bw(base_size = 12) +   ylab("") + xlab('Allen Mouse Cell type') +
  ggtitle('Normalized Tcf4 counts') + 
  theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1,size =8),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()



