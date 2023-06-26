library(tidyverse)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

fn = list.files(pattern = 'idr.optimal_peak.whichtf.tsv', recursive = T, full.names = T)
names(fn) = ss(fn, '/', 2)

alpha = 0.05
whichtf_df = lapply(fn, read_tsv) %>% bind_rows(.id = 'dataset') %>% 
  mutate(logPvalue = Pvalue, 
         Pvalue = 10^-Pvalue, 
         FDR = p.adjust(Pvalue, 'fdr')) %>% 
  filter(FDR < alpha) %>% arrange(Pvalue)

whichtf_df2 = whichtf_df %>% group_by(TF) %>% 
  summarise(n = n(), 
            maxPvalue = max(Pvalue), 
            maxFDR = max(FDR), 
            datasets = paste(dataset, collapse = ', ')) %>% 
  arrange(desc(n), maxPvalue)

output_merged = 'merged.whichtf.xlsx'
writexl::write_xlsx(whichtf_df2, output_merged)

