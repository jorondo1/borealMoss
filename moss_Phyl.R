library(pacman)
p_load(ape, phyloseq, tidyverse, magrittr, ggtree, ggtreeExtra, ggnewscale)
source("myFunctions.R")

# We go through a phyloseq object to easily link trees and MAGs together
# Import PhyloPhlAn output into a (subset) phyloseq object
MAGs.ps <- read.tree("data/RAxML_bestTree.genomes_refined.tre") %>% 
  phy_tree %>% phyloseq %>% 
  #merge with our original phyloseq object; this creates a subset
  merge_phyloseq(readRDS("data/psMossMAGs.RDS")) 

genomes <- MAGs.ps %>% tax_table %>% rownames

# Use the species labels

#Extract tree
tree <- phy_tree(MAGs.ps)

# Create MAG \ species association table
speciesMAG <- MAGs.ps %>% 
  tax_table %>% as.data.frame %>% 
  select(Species) %>% rownames_to_column("MAG")

# Replace tree labels with Species
# MAGs.ps@phy_tree$tip.label <- speciesMAG %$% Species

# using tree as ggtree argument causes FATAL CRASH, using ps for now
p <- ggtree(MAGs.ps, layout="circular") +  # Override the colour mapping shape by creating sham geom_point
  geom_tippoint(mapping = aes(color = Phylum), size = 3.5) 

# Add LFC (from DAA) to significant species
speciesLFC <- readRDS("data/DA_results.RDS") %>% 
  transmute(LFC = lfc_CompartmentGreen, 
            Species = taxon) %>% 
  left_join(speciesMAG, .,
            by = 'Species')
  
hm.mx <- read_tsv("data/novel_quality_scores.txt",
                  col_names = c("MAG", 'comp', 'cont', 'QS')) %>% 
  mutate(MAG = str_remove(MAG, ".fa")) %>% 
  left_join(speciesLFC, by = 'MAG') %>% 
  filter(MAG %in% genomes) %>% 
  column_to_rownames("MAG")

pQS <- gheatmap(p, data = hm.mx["LFC"], 
                offset=0, width=.1, colnames = FALSE,
                ) +
  scale_fill_gradient2(low = 'tan4', high = 'springgreen3', 
                       mid = 'white', na.value = 'grey95',
                      name="Log-fold change \ndifferential abundance",
                       ) + new_scale_fill()

gheatmap(pQS, data = hm.mx["QS"], 
         offset=0.1, width=.1, colnames = FALSE) +
  scale_fill_viridis_c(option="A", name="Quality score")

