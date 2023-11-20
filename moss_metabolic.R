library(pacman); p_load(tidyverse, magrittr, ComplexHeatmap, phyloseq, microViz)
source("moss_functions.R")

# Metabolism
read_tsv("data/microbeannotator_out/metabolic_summary__module_completeness.tab") %>% 
  rename_with(~str_remove_all(.x, "\\.faa\\.ko")) %>% 
  rename_with(~simplify_name(.x))




top_Brown <- df_fig %>% 
  arrange(lfc_CompartmentGreen) %>% 
  slice_head(n = 4) %>% 
  transmute(taxon = taxon,
            lfc = lfc_CompartmentGreen)

top_Green <- df_fig %>% 
  arrange(lfc_CompartmentGreen) %>% 
  slice_tail(n = 4) %>% 
  transmute(taxon = taxon,
            lfc = lfc_CompartmentGreen)






