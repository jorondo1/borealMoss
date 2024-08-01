################################################################################
### AUTHOR : JONATHAN RONDEAU-LECLAIRE #########################################
### Process BLAST output #######################################################
################################################################################

library(pacman)
p_load(tidyverse, magrittr, taxize, doParallel, foreach, RColorBrewer)
options(ENTREZ_KEY = "b69981923f198d72f31791714be0639a2507")
source("scripts/myFunctions.R")

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

### Automatic querying (heavy process)

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

########################################
### Sequence distribution by Kingdom ####
### Fig. S4 data prep

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

# export the essential data (keep this as light as possible)
blastout_short %>% 
  dplyr::select(Host, superkingdom, phylum, Compartment, bp) %>% 
  write_rds("data/R_out/blastout_short.RDS")


########################################################################
### Contamination of P.commune sequences in our MAGs ###################
########################################################################

raw_contam <- read_delim('data/BLAST/all_contam.txt', col_names=c("MAG", "Sample", "MAG_reads", "Host_reads"))

MAG_names <- moss.ps@tax_table %>% as.data.frame %>%  
  dplyr::filter(str_detect(Species,"MAG")) %>% rownames

contamination_summary <- raw_contam %>% 
  group_by(MAG) %>% 
  filter(MAG_reads!=0) %>% 
  filter(MAG %in% MAG_names) %>% 
  summarise(MAG_reads = sum(MAG_reads),
            Host_reads = sum(Host_reads)) %>% 
  mutate(Proportion=100*Host_reads/MAG_reads) %>% 
  arrange(desc(Proportion)) %>% ungroup

contamination_summary

contamination_summary %>%   
  filter(MAG!="Green_N.bin.3") %>% 
  summarise(mean = mean(Proportion),
            sd = sd(Proportion))






