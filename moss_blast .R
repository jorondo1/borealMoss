library(pacman)
p_load(tidyverse, magrittr, taxize, doParallel, foreach, RColorBrewer)
options(ENTREZ_KEY = "b69981923f198d72f31791714be0639a2507")
source("myFunctions.R")

# Parse blast output
colNames<- c("staxids", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", "qseqid", "qlen", "sseqid")
blastout <- Sys.glob("data/blast/*.blastout") %>% 
  map_dfr(~ read.table(.x, sep="\t", col.names = c("staxids", "slen")) %>% 
            mutate(sample = basename(.x)), .x = .) %>% 
  mutate(sample = str_remove(sample, "_subset.blastout"))

# Show top 1000 taxa per sample
topHits <- blastout %>% 
  group_by(staxids, sample) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  ungroup %>% group_by(sample) %>% 
  slice_head(n=1000)

# topBlast <- blastout %>%  filter(staxids %in% topTax)
topTax <- topHits %$% staxids %>% unique

### Automatic querying 
# IF ALREADY DONE : ncbi_tax <- readRDS("data/R_out/taxid_results_full.RDS")

# don't use parallel processing, it generates error because it queries ncbi too much
ncbi_tax <- lapply(
  topTax, # Extract list of unique taxids
  function(id) { tryCatch( { classification(id, db = "ncbi") }, 
                           error = function(e) NA)})
# sometimes querying a lot causes errors, which will introduce NAs in our set.
# Find 
NA_index <- sapply(ncbi_tax, function(x) any(is.na(x))) %>% which
# and correct these :
for (idx in NA_index) {
  ncbi_tax[[idx]] <- classification(topTax[idx], db = 'ncbi')
  # if persistent, remove them :
  if (ncbi_tax[[idx]][1] %>% is.na) {
    ncbi_tax[[idx]] <- NULL
  } # Save that cause it takes time :
} # write_rds(ncbi_tax, "data/R_out/taxid_results_full.RDS")

### Parsing the ncbi query output
taxLvls <- c("staxids","superkingdom", "phylum","class",
             "order","family","genus","species")
# Create an empty df
taxonomy <- data.frame(matrix(ncol = length(taxLvls), nrow = 0))
colnames(taxonomy) <- taxLvls

# Parse the list into the new df:
for (i in 1:length(ncbi_tax)) {
  #df$taxid <- names(results[[i]][1])  # Add taxid as a new column
  taxonomy[i,] <- ncbi_tax[[i]][1][[1]] %>% 
    dplyr::filter(rank %in% taxLvls) %$% name %>% 
    c(names(ncbi_tax[[i]][1]),.) %>% 
  rbind(taxonomy)
}

#################################################
### PLOT 1 : Sequence distribution by Kingdom ####
#################################################

# Create a reduced dataset limited to top species by sample
blastout_short <- blastout %>% 
  # Keep only staxids when they are in the top of a given sample
  semi_join(topHits, by = c("staxids", "sample")) %>% 
  # count number of bases attributed to a taxid, by sample
  group_by(sample, staxids) %>% 
  summarise(bp = sum(slen)) %>% 
  # add Taxonomy
  left_join(taxonomy ,#%>% select(staxids, phylum, superkingdom), 
            by = 'staxids') %>%
  filter(superkingdom %in% c('Bacteria', 'Eukaryota')) %>% 
  # add Sample metadata
  left_join(readRDS('data/R_out/mossMAGs.RDS') %>% 
              sample_data %>% data.frame %>% 
              rownames_to_column('sample') %>% 
              dplyr::select(sample, Host, Location, Compartment),
            by = 'sample') %>% ungroup %>% 
  mutate(Compartment = factor(Compartment, levels = c("Green","Brown")))

blastKingdom <- blastout_short %>% 
  group_by(Host, superkingdom, Compartment) %>% 
  summarise(sum = sum(bp)) %>% 
  slice_head(n=10) %>% 
  ggplot(aes(x = Host, y = sum, fill = superkingdom)) +
  geom_col(position = 'fill') +
  facet_grid('Compartment', scales = 'free') +
  theme_light() + 
  labs(title = 'Proportion of sequence hits',
       y = 'Proportion of BLASTn hits') +
  theme(legend.position = 'bottom')

##################################################
### PLOT 2 : ... by eukaryotic Phylum, by Host ####
##################################################

# Look at what phyla of eukaryotes there are
blast_euk_host <- blastout_short %>% 
  filter(superkingdom == 'Eukaryota') %>% 
  group_by(Host, phylum) %>% 
  summarise(sum = sum(bp)) %>% 
  slice_max(n=8, order_by = sum)

# Custom palette
phylaCols <- colorRampPalette(
  brewer.pal(8, "Set1"))(blast_euk_host$phylum %>% 
                           unique %>% length)
blast_euk_host %>% 
  ggplot(aes(x = Host, y = sum, fill = phylum)) +
  geom_col(position = 'fill') +
  labs(title = 'Eukaryota phyla sequences') +
  scale_fill_manual(values = phylaCols)

########################################################
### PLOT 3 : ...by eukaryotic Phylum and compartment ####
########################################################

blast_euk_sample <- blastout_short %>% 
  filter(superkingdom == 'Eukaryota') %>% 
  group_by(sample, phylum, Compartment, Host) %>% 
  summarise(sum = sum(bp)) %>% 
  ungroup %>% group_by(sample) %>% 
  slice_max(n=4, order_by = sum) 

phylaCols2 <- colorRampPalette(
  brewer.pal(8, "Set1"))(blast_euk_sample$phylum %>% 
                           unique %>% length)
blast_euk_sample %>% 
  ggplot(aes(x = sample, y = sum, fill = phylum)) +
  geom_col(position = 'fill') +
#  facet_wrap('Compartment', ncol = 1, scales = 'free_x') +
  facet_wrap(~ Compartment+Host, nrow=2, scales = 'free_x') +
  labs(title = 'Eukaryota phyla sequences') +
  scale_fill_manual(values = phylaCols2) + 
  theme_light() +
  theme(axis.text.x = element_blank())

########################################################
### PLOT 3 : ...by streptophyta class ####
########################################################

blastout_short %>% 
  filter(phylum == 'Streptophyta') %>% 
  group_by(sample, class, Compartment, Host) %>% 
  summarise(sum = sum(bp)) %>% 
  ungroup %>% group_by(sample) %>% 
  slice_max(n=4, order_by = sum) %>% 
  mutate(Compartment = factor(Compartment, levels = c("Green","Brown"))) %>% 

  ggplot(aes(x = sample, y = sum, fill = class)) +
    geom_col(position = 'fill') +
    facet_wrap(~ Compartment+Host, nrow=2, scales = 'free_x') +
    labs(title = 'Streptophyta class sequences') +
    scale_fill_manual(values = phylaCols2) + 
    theme_light() +
    theme(axis.text.x = element_blank())

