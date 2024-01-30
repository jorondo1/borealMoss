library(pacman)
p_load(tidyverse, magrittr, RColorBrewer, phyloseq, ape)
source("myFunctions.R")

# Sarah's PS object (based on Kraken taxonomic assignment)
sarah.ps<- readRDS("data/ps_comptype.RDS")

###################
#### Abundance #####
###################

abund_GTDB <- parse_SM("data/SM_abund/*gtdb_gather.csv")
abund_MAGs <- parse_SM("data/SM_abund/*custom_gather.csv")
abund_GENB <- parse_SM("data/SM_abund/*genbank_gather.csv")

##################
#### Taxonomy #####
##################
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

# Some taxonomic levels have redundancies because higher levels use alternative names,
# herego there can be a duplicate Order whose Class or Phylum is different. Some 
# can be heterotypic synonyms, others outdated taxonomic names.

# Find duplicates; this function will open the viewer with the names to be resolved.
listDupl <- function(tax, level) {
  level_sym <- rlang::sym(level)
  
  subset <- tax %>% dplyr::select(Domain:!!level_sym) %>% 
    dplyr::group_by(!!level_sym) %>%
    unique
  
  dupList <- subset %>% 
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n > 1) %>%
    dplyr::pull(!!level_sym)
    
  subset %>% arrange(!!level_sym) %>% 
    filter(!!level_sym %in% dupList) %>% print(n = 1000)
}
listDupl(taxonomy, "Genus")

# Correct the different level duplicates, where "New_name" = "Old_name"
corrPhylum <- c("Proteobacteria" = "Pseudomonadota",
               "Actinomycetota" = "Actinobacteriota",
               "Cyanobacteriota" = "Cyanobacteria")
corrClass <- c("Terriglobia" = "Acidobacteriae")
corrOrder <- c("Terriglobales" = "Acidobacteriales",
               "Enterobacterales" = "Enterobacterales_A")
corrFamily <- c("Enterobacteriaceae" = "Enterobacteriaceae_A",
                "SZAS-18" = "UBA10450",
                "Burkholderiaceae" = "Burkholderiaceae_B")
# correct, but check family too...
taxonomy %<>% mutate(Phylum = recode(Phylum, !!!corrPhylum),
                    Class = recode(Class, !!!corrClass),
                    Order = recode(Order, !!!corrOrder),
                    Family = recode(Family, !!!corrFamily))

##################
#### Metadata #####
##################

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
psMossMAGs <- makePhyloSeq(abund_MAGs, sampleData, taxonomy, 
                           tree = "data/RAxML_bestTree.genomes_refined.tre") %>% 
  # TEMPORARY : REMOVE MAG identified as contaminated by GUNC
  prune_taxa(taxa_names(.) !="Green_AO.bin.6",.)

write_rds(psMossGTDB,"data/R_out/psMossGTDB.RDS")
write_rds(psMossMAGs,"data/R_out/psMossMAGs.RDS")

