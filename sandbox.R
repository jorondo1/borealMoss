library(pacman)
p_load(dplyr, vegan, phyloseq )
moss.ps <- readRDS("data/R_out/mossMAGs.RDS")

# Library size
read_delim("data/raw_sequence_counts.txt", col_names = F) %>% 
  summarise(mean_seq = mean(X3),
            sd_seq = sd(X3))

read_delim("data/sequence_counts.txt", col_names = F) %>%
  summarise(mean_seq = mean(X2),
            sd_seq = sd(X2))

# Check eukaryota proportion using Genbank data
GBNK_melt <- mossGBNK.ps %>% psmelt

GBNK_melt %>% group_by(Sample, Domain, Compartment, Host) %>% 
  relab.fun %>% 
  group_by(Host, Compartment, Domain) %>% 
  summarise(meanRelativeAbund = mean(relAb))

# Ok they are not negligible. who are they ?
GBNK_melt %>% group_by(Sample, Domain, Phylum, Host) %>% 
  relab.fun %>% 
  group_by(Host, Domain, Phylum) %>% 
  summarise(meanRelativeAbund = mean(relAb),
            SD = sd(relAb)) %>% 
  filter(meanRelativeAbund > 0.01 & Domain == 'Eukaryota')

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
