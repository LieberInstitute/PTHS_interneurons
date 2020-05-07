library(tidyverse)
library(ggplot2)
library(SingleCellExperiment)
library(DESeq2)
library(ggplot2)

options(stringsAsFactors = F)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

enrichFisher = function(genes, markerList, universe){
  genes = genes[genes %in% universe]
  fisherMat = do.call('rbind',lapply(markerList, function(markers){
    markers = markers[markers %in% universe]
    test = fisher.test(table(universe %in% markers, universe %in% genes))
    return(data.frame(p.value = test$p.value,
                      OR = test$estimate,
                      OR.conf.int.max = max(test$conf.int),
                      OR.conf.int.min = min(test$conf.int)))
  }))
  return(fisherMat)
}

##################################################
# load in Tasic cell type specific marker genes #
load('rdas/Tasic2018_allen_mouse_markerGenes.rda')

# filter markers by FDR < 0.05
alpha = 0.05
markers.All = lapply(markers.tTestUpAll, function(x) as.data.frame(x) %>% rownames_to_column(var = "Gene") %>% filter(FDR < alpha))
markers.Some = lapply(markers.tTestUpSome, function(x) as.data.frame(x) %>% rownames_to_column(var = "Gene") %>% filter(FDR < alpha))

# drop celltypes w/ 0 marker genes
indKeep.AllMarkers = which(sapply(markers.All, nrow) > 0)
indKeep.SomeMarkers = which(sapply(markers.Some, nrow) > 0)
markerList.All = lapply(markers.All[indKeep.AllMarkers],'[[', 'Gene')
markerList.Some = lapply(markers.Some[indKeep.SomeMarkers],'[[', 'Gene')
names(markerList.All) = paste('UpInAll',names(markerList.All),sep = '.')
names(markerList.Some) = paste('UpInSome',names(markerList.Some),sep = '.')
markerList = c(markerList.All, markerList.Some)


###################################
# load in PTHS mouse mega DESeq object #
load('/users/bphan/tcf4/PTHS_mouse/tcf4_mouse/rdas/mega_tcf4_ages_DE_objects_DESeq2.rda', 
     envir = mega <- new.env())

# background genes
universe = intersect(rownames(markers.tTestUpAll[[1]]), mega$outGeneList[[1]]$Symbol)
length(universe) #20282 genes in both datasets

# filter Phan et al. DEGs
mega$sigGeneList2 = c(lapply(mega$sigGeneList, function(x) x %>% filter(!is.na(padj), padj < alpha, log2FoldChange > 0)),
                 lapply(mega$sigGeneList, function(x) x %>% filter(!is.na(padj), padj < alpha, log2FoldChange < 0)))
# mega$sigGeneList2 = c(lapply(mega$sigGeneList, function(x) x %>% filter(log2FoldChange > 0)),
#                  lapply(mega$sigGeneList, function(x) x %>% filter(log2FoldChange < 0)))
names(mega$sigGeneList2) = paste(names(mega$sigGeneList2), rep(c('Up','Down'), each = length(mega$sigGeneList)),sep = '.')
degList = lapply(mega$sigGeneList2,'[[','Symbol')
names(degList) = paste('mega',names(degList), sep = '.')

# perform Fisher Overlap enrichment test w/ All markers
enrichMat1 = do.call('rbind',lapply(degList, enrichFisher, markerList = markerList.All, universe = universe))
enrichMat1$pthsGroup = ss(rownames(enrichMat1),'\\.', 1)
enrichMat1$ageGroup = ss(rownames(enrichMat1),'\\.', 2)
enrichMat1$DEGdirGroup = factor(ss(rownames(enrichMat1),'\\.', 3),levels = c('Up','Down'))
enrichMat1$cutoffGroup = ss(rownames(enrichMat1),'\\.', 4)
enrichMat1$markerGroup = ss(rownames(enrichMat1),'\\.', 5)
enrichMat1 = enrichMat1 %>% mutate(FDR = p.adjust(p.value, 'BH'))
signifEnrichMat1 = enrichMat1 %>% filter(FDR < alpha & OR > 1)
# table(signifEnrichMat1$ageGroup, signifEnrichMat1$markerGroup, signifEnrichMat1$cutoffGroup)

#######################
## reorder celltypes ##
keepColNames = c('class_label','class_color','class_order', 'subclass_label','subclass_color','subclass_order')
cellData = read.csv('data/sample_annotations.csv') %>%  select(all_of(keepColNames)) %>% distinct() %>% 
  filter(class_label !='Exclude') %>% arrange(class_order, subclass_order)
cellData$subclass_label = factor(cellData$subclass_label, levels = cellData$subclass_label)
rownames(cellData) = cellData$subclass_label

###########################
## make enrichment plots ##
dat = cbind(enrichMat1, cellData[enrichMat1$markerGroup,])
dat$FDR_signif = dat$FDR < alpha
dat$logOR = log10(dat$OR)
dat$markerGroup = factor(dat$markerGroup, levels = rev(cellData$subclass_label))

subclass_color = cellData$subclass_color
names(subclass_color) = cellData$subclass_label

pdf('plots/Tasic_2018_Allen_mouse_markerEnrichment_PTHSmouse.pdf', height = 8, width = 5)
ggplot(subset(dat, cutoffGroup =='UpInAll' & ageGroup =='Adult'), 
       aes(x = markerGroup, y = OR, fill = subclass_label, alpha = FDR_signif, color = FDR_signif)) +
  geom_bar(stat = 'identity') + coord_flip() + scale_y_log10() +
  geom_errorbar(aes(ymin = OR.conf.int.min, ymax = OR.conf.int.max), width=0.2) +
  facet_grid(class_label ~ DEGdirGroup, scales = 'free_y', space = 'free') + 
  scale_fill_manual( values = subclass_color, guide = FALSE) + 
  scale_alpha_manual( values = c(0.5,1), guide = FALSE) + 
  scale_color_manual( values = c(NA,'black'), guide = FALSE) + 
  geom_hline(yintercept = 1, color = 'black') + 
  xlab('Cell type') + ylab('Odds Ratio') + 
  theme_bw(base_size = 12) + 
  ggtitle('Adult DEGs')


ggplot(subset(dat, cutoffGroup =='UpInAll' & ageGroup =='p1'), 
       aes(x = markerGroup, y = OR, fill = subclass_label, alpha = FDR_signif, color = FDR_signif)) +
  geom_bar(stat = 'identity') + coord_flip() + scale_y_log10() +
  geom_errorbar(aes(ymin = OR.conf.int.min, ymax = OR.conf.int.max), width=0.2) +
  facet_grid(class_label ~ DEGdirGroup, scales = 'free_y', space = 'free') + 
  scale_fill_manual( values = subclass_color, guide = FALSE) + 
  scale_alpha_manual( values = c(0.5,1), guide = FALSE) + 
  scale_color_manual( values = c(NA,'black'), guide = FALSE) + 
  geom_hline(yintercept = 1, color = 'black') + 
  xlab('Cell type') + ylab('Odds Ratio') + 
  theme_bw(base_size = 12) + 
  ggtitle('P1 DEGs')
dev.off()



###################################
# load in PTHS mouse by model DESeq object #
load('/users/bphan/tcf4/PTHS_mouse/tcf4_mouse/rdas/mega_tcf4_separate_DE_objects_DESeq2.rda', 
     envir = separate <- new.env())

# filter Phan et al. DEGs of each model
separate$sigGeneList2 = c(lapply(separate$sigGeneList, function(x) x %>% filter(!is.na(padj), padj < alpha, log2FoldChange > 0)),
                      lapply(separate$sigGeneList, function(x) x %>% filter(!is.na(padj), padj < alpha, log2FoldChange < 0)))
names(separate$sigGeneList2) = paste(names(separate$sigGeneList2), rep(c('Up','Down'), each = length(mega$sigGeneList)),sep = '.')
degList2 = c(degList, lapply(separate$sigGeneList2,'[[','Symbol'))
degList2 = degList2[grepl('mega|Maher', names(degList2))]
degList2 = degList2[grepl('Adult', names(degList2))]

# enrichment 
enrichMat2 = do.call('rbind',lapply(degList2, enrichFisher, markerList = markerList.All, universe = universe))
enrichMat2$pthsGroup = factor(ss(rownames(enrichMat2),'\\.', 1), 
                              levels = c('mega','Maher'))
enrichMat2$ageGroup = ss(rownames(enrichMat2),'\\.', 2)
enrichMat2$DEGdirGroup = factor(ss(rownames(enrichMat2),'\\.', 3),levels = c('Up','Down'))
enrichMat2$cutoffGroup = ss(rownames(enrichMat2),'\\.', 4)
enrichMat2$markerGroup = ss(rownames(enrichMat2),'\\.', 5)
enrichMat2 = enrichMat2 %>% mutate(FDR = p.adjust(p.value, 'BH'))

# plot enrichments 
dat2 = cbind(enrichMat2, cellData[enrichMat2$markerGroup,])
dat2$FDR_signif = dat2$FDR < alpha
dat2$logOR = log10(dat2$OR)
dat2$markerGroup = factor(dat2$markerGroup, levels = rev(cellData$subclass_label))

pdf('plots/Tasic_2018_Allen_mouse_markerEnrichment_PTHSmouse_separate_models.pdf', width = 8, height = 8)
ggplot(subset(dat2, cutoffGroup =='UpInAll' & ageGroup =='Adult'), 
       aes(x = markerGroup, y = OR, fill = subclass_label, alpha = FDR_signif, color = FDR_signif)) +
  geom_bar(stat = 'identity') + coord_flip() + scale_y_log10() +
  geom_errorbar(aes(ymin = OR.conf.int.min, ymax = OR.conf.int.max), width=0.2) +
  facet_grid(class_label ~ DEGdirGroup + pthsGroup, scales = 'free_y', space = 'free') + 
  scale_fill_manual( values = subclass_color, guide = FALSE) + 
  scale_alpha_manual( values = c(0.5,1), guide = FALSE) + 
  scale_color_manual( values = c(NA,'black'), guide = FALSE) + 
  geom_hline(yintercept = 1, color = 'black') + 
  xlab('Cell type') + ylab('Odds Ratio') + 
  theme_bw(base_size = 12) + 
  ggtitle('Adult DEGs')


ggplot(subset(dat2, cutoffGroup =='UpInAll' & ageGroup =='Adult'), 
       aes(x = DEGdirGroup, y = OR, fill = pthsGroup, alpha = FDR_signif, color = FDR_signif)) +
  geom_bar(stat = 'identity', position = 'dodge') + coord_flip() + scale_y_log10() +
  geom_errorbar(aes(ymin = OR.conf.int.min, ymax = OR.conf.int.max), width=0.2, 
                position = position_dodge(width = .9)) +
  facet_wrap(~markerGroup, scales = 'free') + 
  #scale_fill_manual( values = subclass_color, guide = FALSE) + 
  scale_alpha_manual( values = c(0.5,1), guide = FALSE) + 
  scale_fill_brewer(palette = 'Set1') + 
  scale_color_manual( values = c(NA,'black'), guide = FALSE) + 
  geom_hline(yintercept = 1, color = 'black') + 
  xlab('Cell type') + ylab('Odds Ratio') + 
  theme_bw(base_size = 12) + 
  ggtitle('Adult DEGs')
dev.off()


ggplot(subset(dat2, cutoffGroup =='UpInAll' & ageGroup =='p1'), 
       aes(x = markerGroup, y = OR, fill = subclass_label, alpha = FDR_signif, color = FDR_signif)) +
  geom_bar(stat = 'identity') + coord_flip() + scale_y_log10() +
  geom_errorbar(aes(ymin = OR.conf.int.min, ymax = OR.conf.int.max), width=0.2) +
  facet_grid(class_label ~ DEGdirGroup, scales = 'free_y', space = 'free') + 
  scale_fill_manual( values = subclass_color, guide = FALSE) + 
  scale_alpha_manual( values = c(0.5,1), guide = FALSE) + 
  scale_color_manual( values = c(NA,'black'), guide = FALSE) + 
  geom_hline(yintercept = 1, color = 'black') + 
  xlab('Cell type') + ylab('Odds Ratio') + 
  theme_bw(base_size = 12) + 
  ggtitle('P1 DEGs')
dev.off()
