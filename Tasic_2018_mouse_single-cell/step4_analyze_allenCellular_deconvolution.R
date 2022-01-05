library(tidyverse)
library(DESeq2)
library(reshape2)
library (MGLM)
library(ggplot2)

options(stringsAsFactors = FALSE)

#######################
## reorder celltypes ##
keepColNames = c('class_label','class_color','class_order', 'subclass_label','subclass_color','subclass_order')
cellData = read.csv('data/sample_annotations.csv') %>%  select(all_of(keepColNames)) %>% distinct() %>% 
  filter(class_label !='Exclude') %>% arrange(class_order, subclass_order)
cellData$subclass_label = make.names(cellData$subclass_label)
cellData$subclass_label = factor(cellData$subclass_label, levels = cellData$subclass_label)
rownames(cellData) = cellData$subclass_label

subclass_color = cellData$subclass_color
names(subclass_color) = cellData$subclass_label

###################################
# load in PTHS mouse mega DESeq object #
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mega_tcf4_ages_DESeq2_svaAdj.rda')
load('/dcl01/lieber/ajaffe/Brady/mouseRNAseq/mega_tcf4_ages_DESeq2_svaAdj.rda')

# get phenotype data
shared_columns = intersect(names(colData(geneDds$p1)), names(colData(geneDds$Adult)))
pd =as.data.frame( rbind(colData(geneDds$p1)[,shared_columns], 
                         colData(geneDds$Adult)[,shared_columns]))
pd$Line = factor(pd$Line, levels =  c('Maher','Act','Nest','Del','R579W','Sweatt'))

# get data.frame of deconvolved cell proportions
df = read.csv("tables/PTHSmouse_logCounts_deconvolution_score_p1-Adult_maxIter10000.csv", row.names = 1)
df = df[, apply(df,2, mean) > 0.001] # drop cell types with low predicted cell type proportions
df = df/ rowSums(df) # renormalize proportions
all.equal(rownames(df), rownames(pd))

dat = cbind(pd, df)
datLong = melt(dat, id.vars = names(pd), value.name = 'Proportion' )
datLong$class_label = cellData[as.character(datLong$variable), 'class_label']
datLong$variable = factor(datLong$variable, levels = rev(cellData$subclass_label))


pdf('plots/PTHSmouse_logCounts_deconvolution_GABAergic.pdf', width = 9, height = 11)
ggplot(subset(datLong, class_label=='GABAergic' & Age == 'Adult'), 
       aes(x = Line, y = Proportion, fill = Genotype)) + 
  geom_boxplot() + geom_point(position = position_jitterdodge(), pch = 21) + 
  facet_wrap(~variable, scales = 'free_y', ncol = 2) + 
  scale_fill_manual(values = c('darkgray','red')) + 
  theme_bw(base_size = 14)
dev.off()

pdf('plots/PTHSmouse_logCounts_deconvolution_Glutamatergic.pdf',width = 9, height= 20)
ggplot(subset(datLong, class_label=='Glutamatergic'& Age == 'Adult'), 
       aes(x = Line, y = Proportion, fill = Genotype)) + 
  geom_boxplot() + geom_point(position = position_jitterdodge(), pch = 21) + 
  facet_wrap(~variable, scales = 'free_y', ncol = 3) + 
  scale_fill_manual(values = c('darkgray','red')) + 
  theme_bw(base_size = 14)
dev.off()

pdf('plots/PTHSmouse_logCounts_deconvolution_Non-neuronal.pdf',width = 9, height= 5)
ggplot(subset(datLong, class_label=='Non-neuronal'& Age == 'Adult'), 
       aes(x = Line, y = Proportion, fill = Genotype)) + 
  geom_boxplot() + geom_point(position = position_jitterdodge(), pch = 21) + 
  facet_wrap(~variable, scales = 'free_y', ncol = 3) + 
  scale_fill_manual(values = c('darkgray','red')) + 
  theme_bw(base_size = 14)
dev.off()



######################
# only look at Maher Tcf4-tr line 
# get data.frame of deconvolved cell proportions
indKeep = which(pd$Line =='Maher')
df2 = df[rownames(pd)[indKeep], apply(df,2, mean) > 0.005] # drop cell types with low predicted cell type proportions
df2 = df2/ rowSums(df2) # renormalize proportions

dat2 = cbind(pd[indKeep,], df2)
dat2Long = melt(dat2, id.vars = names(pd), value.name = 'Proportion' )
dat2Long$class_label = cellData[as.character(dat2Long$variable), 'class_label']
dat2Long$variable = factor(dat2Long$variable, levels = rev(cellData$subclass_label))

pdf('plots/PTHSmouse_logCounts_deconvolution_GABAergic_Tcf4-trLine.pdf',width = 12, height= 6)
ggplot(subset(dat2Long, Line=='Maher' & class_label =='GABAergic' ), 
       aes(x = Age, y = Proportion, fill = Genotype)) + 
  geom_boxplot() + geom_point(position = position_jitterdodge(), pch = 21) + 
  facet_wrap(~variable, scales = 'free', ncol = 4) + 
  scale_fill_manual(values = c('darkgray','red')) + 
  theme_bw(base_size = 14) + 
  xlab('Cell subclass')
dev.off()

pdf('plots/PTHSmouse_logCounts_deconvolution_Glutamatergic_Tcf4-trLine.pdf',width = 12, height= 9)
ggplot(subset(dat2Long, Line=='Maher' & class_label =='Glutamatergic' ), 
       aes(x = Age, y = Proportion, fill = Genotype)) + 
  geom_boxplot() + geom_point(position = position_jitterdodge(), pch = 21) + 
  facet_wrap(~variable, scales = 'free', ncol = 4) + 
  scale_fill_manual(values = c('darkgray','red')) + 
  theme_bw(base_size = 14) + 
  xlab('Cell subclass')
dev.off()

pdf('plots/PTHSmouse_logCounts_deconvolution_Non-neuronal_Tcf4-trLine.pdf',width = 12, height= 3)
ggplot(subset(dat2Long, Line=='Maher' & class_label =='Non-neuronal' ), 
       aes(x = Age, y = Proportion, fill = Genotype)) + 
  geom_boxplot() + geom_point(position = position_jitterdodge(), pch = 21) + 
  facet_wrap(~variable, scales = 'free', ncol = 4) + 
  scale_fill_manual(values = c('darkgray','red')) + 
  theme_bw(base_size = 14) + 
  xlab('Cell subclass')
dev.off()


# multinomial regression
Y = as.matrix(round(df/min(df[df>0])) + 1)
mnReg <- MGLMreg(Y~ Genotype+Age + Line, pd, dist="DM")
show(mnReg)
# Hypothesis test: 
#   wald value     Pr(>wald)
# (Intercept)  2047.2587  0.000000e+00
# GenotypeHT     88.0024  2.306664e-07
# AgeAdult     2083.6182  0.000000e+00
# LineAct       410.9342  6.431604e-68
# LineNest      646.7405 2.803867e-116
# LineDel       774.0838 8.395224e-143
# LineR579W     668.2672 9.524199e-121
# LineSweatt    678.6203 6.717462e-123


cellTypes = colnames(df)
names(cellTypes) =colnames(df)

########################################################
# test change in proportions of cell types across all ages
lmList = lapply(cellTypes, function(i){
  reg1 = lm(df[,i] ~ Genotype  + Line + Age, data = pd)
  summary(reg1)
})
lmList[['Pvalb']]
p.value =  sapply(lmList, function(x) x$coefficients['GenotypeHT','Pr(>|t|)'])
FDR = p.adjust(p.value,'BH') 
p.value[order(p.value)]
FDR[order(FDR)]


########################################################
# test change in proportions of cell types in Adults only
lmList2 = lapply(cellTypes, function(i){
  indKeep = which(pd$Age =='Adult')
  reg2= lm(df[indKeep,i] ~ Genotype + Line , data = pd[indKeep,])
  summary(reg2)
})
lmList2[['Pvalb']]
p.value2 =  sapply(lmList2, function(x) x$coefficients['GenotypeHT','Pr(>|t|)'])
FDR2 = p.adjust(p.value2,'BH')
FDR2[order(FDR2)]



########################################################
# test change in proportions of cell types in Adults only
lmList3 = lapply(cellTypes, function(i){
  indKeep = which(pd$Age =='p1')
  reg2= lm(df[indKeep,i] ~ Genotype + Line , data = pd[indKeep,])
  summary(reg2)
})
lmList3[['Pvalb']]
p.value3 =  sapply(lmList3, function(x) x$coefficients['GenotypeHT','Pr(>|t|)'])
FDR3 = p.adjust(p.value3,'BH')
FDR3[order(FDR3)]


pdf('plots/PTHSmouse_logCounts_deconvolution_Oligo_Pvalb.pdf', width = 7, height = 5)
ggplot(subset(datLong, variable %in% c('Pvalb','Oligo')), 
       aes(x = Line, y = Proportion, fill = Genotype)) + 
  geom_boxplot() + geom_point(position = position_jitterdodge(), pch = 21) + 
  facet_grid(variable~Age, scales = 'free', space = 'free') + 
  scale_fill_manual(values = c('darkgray','red')) + 
  theme_bw(base_size = 12)
dev.off()

pdf('plots/PTHSmouse_logCounts_deconvolution_Pvalb_stacked.pdf', width = 3.5, height = 5)
ggplot(subset(datLong, variable %in% c('Pvalb')), 
       aes(x = Line, y = Proportion * 100, fill = Genotype)) + 
  geom_boxplot() + geom_point(position = position_jitterdodge(), pch = 21) + 
  facet_wrap(~Age, scales = 'free', ncol= 1) + 
  scale_fill_manual(values = c('darkgray','red')) + 
  theme_bw(base_size = 12) + 
  xlab('Mouse Line') + ylab('Est. % of Pvalb+ Interneurons') + 
  theme(legend.position = 'bottom')
dev.off()

