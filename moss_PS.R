library(pacman)
p_load(tidyverse, magrittr, RColorBrewer, phyloseq, ape)
source("moss_functions.R")

# Sarah's PS object (based on Kraken taxonomic assignment)
sarah.ps<- readRDS("ps_comptype.RDS")
abund_GTDB <- parse_SM("data/SM_abund/*gtdb_gather.csv")
abund_MAGs <- parse_SM("data/SM_abund/*custom_gather.csv")

# GTDB Taxonomy
taxLvls <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

taxonomy <- read_delim("data/gtdbtk_summary.tsv") %>% 
  separate(classification, into = taxLvls, sep = ";", fill = "right") %>% 
  
  # We rename the bins, because we'll add these to the species
  # (since MAG species is not always known, this prevents two species-level bins
  # from the same genus to have the same species name)
  mutate(unique_MAGs = simplify_name(user_genome),
         genome = user_genome) %>% 
  dplyr::select(all_of(c('genome', taxLvls, 'unique_MAGs'))) %>% 
  
  # Adding known GTDB taxa
  bind_rows(
    read_delim("data/gtdb_taxonomy_subset.csv", delim=',', col_names = F) %>% 
      select(-X2) %>% 
      set_names(c('genome',taxLvls))
    ) %>% 
  mutate_at(taxLvls, rename.fun) %>% 
  
  # If the species is unknown, we'll fill the Species field with the nearest known
  # taxonomic rank, and add a unique MAG identifier
  mutate(across(everything(), ~na_if(.x, "")),
    Species = ifelse(
      is.na(Species),
      paste0(coalesce(Genus, Family, Order, Class, Phylum, Domain), 
             " ", unique_MAGs),
      Species
    )) %>% select(-unique_MAGs)

# Species was used to designate Host Species, risk of confusion with Microbiome species
sampleData <- sarah.ps %>% sample_data %>% as("data.frame") %>% 
  mutate(Host = str_c( 
    str_extract(Species, "^[A-Z]"), # Extract the first capital letter
    "_", # add underscore
    str_extract(Species, "(?<= )[a-z]+"), # Extract the second word (species name)
    sep = ""
  ), .keep = 'unused')

# Generate PS objects
psMossGTDB <- makePhyloSeq(abund_GTDB, sampleData, taxonomy)
psMossMAGs <- makePhyloSeq(abund_MAGs, sampleData, taxonomy)

write_rds(psMossGTDB,"data/psMossGTDB.RDS")
write_rds(psMossMAGs,"data/psMossMAGs.RDS")

