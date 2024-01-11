library(pacman)
p_load(ape, phyloseq, tidyverse, magrittr, ggtree, ggtreeExtra, ggnewscale)
source("myFunctions.R")

# We go through a phyloseq object to easily link trees and MAGs together
# Import PhyloPhlAn output into a (subset) phyloseq object
MAGs.ps <- read.tree("data/RAxML_bestTree.genomes_refined.tre") %>% 
  phy_tree %>% phyloseq %>% 
  # merge with our original phyloseq object; this creates a subset
  merge_phyloseq(readRDS("data/psMossMAGs.RDS")) 

#####################################
#### PLOT 1. MAGs characteristics ####
#####################################

# Create MAG \ species association table
speciesMAG <- MAGs.ps %>% 
  tax_table %>% as.data.frame %>% 
  select(Species) %>% rownames_to_column("MAG")

# Replace tree labels with Species
# MAGs.ps@phy_tree$tip.label <- speciesMAG %$% Species

# Add LFC (from DAA) to significant species
speciesLFC <- readRDS("data/DA_results.RDS") %>% 
  transmute(LFC = lfc_CompartmentGreen, 
            Species = taxon) %>% 
  left_join(speciesMAG, .,
            by = 'Species')

# DF for MAG quality heatmap 
MAG_names <- MAGs.ps %>% tax_table %>% rownames # MODIFY LATER to only get names with .bin. 
hm.mx <- read_tsv("data/novel_quality_scores.txt",
                  col_names = c("MAG", 'comp', 'cont', 'QS')) %>% 
  mutate(MAG = str_remove(MAG, ".fa")) %>% 
  left_join(speciesLFC, by = 'MAG') %>% 
  filter(MAG %in% MAG_names) %>% 
  column_to_rownames("MAG")

# using tree as ggtree argument causes FATAL CRASH, using ps for now

# Expand a Brewer palette to 12 colours:
nTax <- MAGs.ps@tax_table %>% as.data.frame %$% Class %>% unique %>% length
myColors <- colorRampPalette(brewer.pal(9, "Set1"))(nTax+1) # +1 for NAs

# Plot the tree
p <- ggtree(MAGs.ps, layout="circular") +  # Override the colour mapping shape by creating sham geom_point
  geom_tippoint(mapping = aes(color = Class), size = 3.5) +
  scale_colour_manual(values = myColors, na.value = "grey10")

# Add a heatmap
gheatmap(p, data = hm.mx["QS"], 
         offset=0.01, width=.1, colnames = FALSE) +
  scale_fill_viridis_c(option="A", name="MAG Quality Score") +
  labs(color = "GTDB-Tk Class Assignment")

#####################################
#### PLOT 1. MAGs characteristics ####
#####################################

# Subset taxa for tree layer
DA_results <- readRDS('data/DA_results.RDS') %>% 
  transmute(Species = taxon,
            compAss = case_when(lfc_CompartmentGreen>0 ~ "Green",
                                    lfc_CompartmentGreen<0 ~ "Brown"))

DA_species.ps <- subset_taxa(MAGs.ps, Species %in% (DA_results %$% Species))

# create table w/ normalised abundance 
relab.df <- MAGs.ps %>% psmelt %>%
  group_by(Sample) %>% 
  mutate(total = sum(Abundance)) %>% 
  transmute(OTU = OTU,
            relAb = (Abundance/total),
            Compart = Compartment,
            Species = Species) %>% 
  inner_join(DA_results, by = "Species") %>%  # add compartment association
  #filter(Compartment == compAss)
  group_by(Compart, OTU, Species, compAss) %>% 
  summarise(mean_relAb = mean(relAb)) %>% arrange(OTU) %>% 
  mutate(across(where(is.character),as.factor)) %>% 
  ungroup
  
# Add a fake taxonomic level to access compAss by species
tax_table(DA_species.ps) <- DA_species.ps@tax_table %>% 
  as.data.frame %>% rownames_to_column("OTU") %>% 
  left_join((relab.df %>% 
              dplyr::select(Species, compAss) %>% 
              unique), by = "Species") %>% 
  column_to_rownames("OTU") %>% as.matrix

ggtree(DA_species.ps,layout = "circular") +
  geom_tippoint(mapping = aes(color = compAss), size = 3.5) +
  scale_colour_manual(values = c('lightsalmon4', 'darkolivegreen')) +
  geom_fruit(relab.df, 
                geom = geom_col, 
                mapping = aes(y = OTU, 
                              x = mean_relAb, 
                              fill = Compart),
                position = position_dodgex(width = 0.7),
                pwidth = .6,
                offset = .2, 
                # Add a grid behind
                axis.params=list(
                  axis       = "x",
                  text.size  = 1.8,
                  hjust      = 1,
                  vjust      = 0.5,
                  nbreak     = 3,
                ),
                grid.params=list()
                ) +
  scale_fill_manual(values = c('lightsalmon4', 'darkolivegreen')) +
  labs(fill = "Sample type relative abundance",
       colour = "Species compartment association") +
  theme(legend.position = "bottom") +
  guides(
    color = guide_legend(title.position = "top"),
    fill = guide_legend(title.position = "top")
  ) 


# plot circular tree
# pQS <- gheatmap(p, data = hm.mx["LFC"], 
#                 offset=0, width=.1, colnames = FALSE,
#                 ) +
# scale_fill_gradient2(low = 'tan4', high = 'springgreen3', 
#                      mid = 'white', na.value = 'grey95',
#                     name="Log-fold change \ndifferential abundance",
#                      ) + new_scale_fill() 

# research credits report figure 
gheatmap(p, data = hm.mx["QS"], 
         offset=0.01, width=.1, colnames = FALSE) +
  scale_fill_viridis_c(option="A", name="Quality score") +
 # ggtitle("Novel species MAGs phylogeny and assembly quality scores.") +
  theme(text = element_text(size = 26))

