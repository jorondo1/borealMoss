library(pacman)
p_load(phyloseq, tidyverse, magrittr)
source("myFunctions.R")

# Merging two MicrobeAnnotator outputs
pw <- full_join(
  read_delim("data/metabolic_summary__module_completeness_missing.tab"),
  read_delim("data/metabolic_summary__module_completeness.tab"),
  by = c('module', 'name', 'pathway group')
  ) %>% 
  dplyr::select(-name, -`pathway group`) %>% 
  rename_with(~str_remove_all(.x, "\\.faa\\.ko")) %>% 
  rename_with(~simplify_name(.x)) # shortens MAG names

# transpose to get MAG name as rows, modules as columns :
pw_t <- pw %>% dplyr::select(-module) %>% t %>% 
  data.frame %>% setNames(pw$module)




# DA results only list the Species, not the shortened MAG name
# Pathways matrix only has shortened MAG name as rows + modules as columns
# I extracted the Species for these MAGs from the taxonomy list in JO_psMossMAGs.RDS file
# But the MAG names written in the taxonomy file are written in the long form, so I had to shorten the names with one of the functions you wrote
# I think I just x_joined (can't remember if it was left or inner) to match the Species to the pathway matrix names then left_joined that to DA results to get a matrix of all the DA output + modules
