library(pacman)
p_load(dplyr, vegan, phyloseq )
moss.ps <- readRDS("data/R_out/mossMAGs.RDS")

# Which species are present in every ....
find_core <- function(ps) {
  IDs <- otu_table(ps) %>%
    vegan::decostand(method = 'total') %>%
    data.frame %>%
    filter(apply(., 1, function(x) all(x != 0))) %>% 
    rownames
  tax_table(ps)[IDs,]
}
subset_samples(moss.ps, Compartment == 'Green') %>% find_core(.)
subset_samples(moss.ps, Compartment == 'Brown') %>% find_core(.)

subset_samples(moss.ps, Host == 'P_commune') %>% find_core(.)



# Lichenibacterium?
moss.ps@tax_table %>% as.data.frame %>% 
  filter(Family == 'Lichenibacterium') #nope
