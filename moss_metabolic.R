library(pacman)
p_load(phyloseq, tidyverse, magrittr)
source("myFunctions.R")

pw <- read_delim("data/microbeannotator_out/metabolic_summary__module_completeness.tab")

done <- pw %>% select(-module, -`pathway group`, -name) %>% 
  colnames %>% str_remove(".faa.ko")

moss.ps %>% prune_taxa(done, .) %>% tax_table %>% View
