library(pacman)
p_load(ape, phyloseq, tidyverse, magrittr, ggtree,
       RColorBrewer, ggtreeExtra, ggnewscale)
source("myFunctions.R")

# We go through a phyloseq object to easily link trees and MAGs together
# Import PhyloPhlAn output into a (subset) phyloseq object
moss.ps <- read.tree("data/RAxML_bestTree.genomes_refined.tre") %>% 
   phy_tree %>% phyloseq %>% 
  # merge with our original phyloseq object; this creates a subset
  merge_phyloseq(readRDS("data/psMossMAGs.RDS")) 

moss.melt <- psmelt(moss.ps)

# Add LFC (from DAA) to significant species
speciesLFC <- readRDS("data/DA_results.RDS") %>% 
  transmute(LFC = lfc_CompartmentGreen, 
            Species = taxon, 
            compAss = case_when(lfc_CompartmentGreen>0 ~ "Green",
                                lfc_CompartmentGreen<0 ~ "Brown")) %>% 
  right_join(moss.ps %>% # identifier \ species association table
              tax_table %>% as.data.frame %>% 
              select(Species) %>% rownames_to_column("MAG"),
            by = 'Species')

#####################################
#### PLOT 1. MAGs characteristics ####
#####################################

# Replace tree labels with Species
# moss.ps@phy_tree$tip.label <- speciesMAG %$% Species

# DF for MAG quality heatmap 
MAG_names <- moss.ps@tax_table %>% rownames %>% .[grep(".bin.", .)]
hm.mx <- read_tsv("data/novel_quality_scores.txt",
                  col_names = c("MAG", 'comp', 'cont', 'QS')) %>% 
  mutate(MAG = str_remove(MAG, ".fa")) %>% 
  left_join(speciesLFC, by = 'MAG') %>% 
  filter(MAG %in% MAG_names) %>% 
  column_to_rownames("MAG")

# Expand a Brewer palette to 12 colours:
nClass <- moss.ps@tax_table %>% as.data.frame %$% Class %>% unique %>% length
classCol <- colorRampPalette(brewer.pal(9, "Set1"))(nClass+1) # +1 for NAs

# Plot the tree
p <- moss.ps %>% prune_taxa(taxa = MAG_names, .) %>% 
  ggtree(layout="fan", size=0.2) +  # Override the colour mapping shape by creating sham geom_point
  xlim(-0.6, NA) + # prevent the high-level branches from clustering in the middle
  geom_tippoint(mapping = aes(color = Class), size = 2.5) +
  scale_colour_manual(values = classCol, na.value = "grey10")

# Add a heatmap
###! !!!!! Check if QS is aligned, there's no anchoring!
gheatmap(p, data = hm.mx["QS"], 
         offset=0.01, width=.1, colnames = FALSE) +
  scale_fill_viridis_c(option="A", name="MAG Quality Score") +
  labs(color = "GTDB-Tk Class Assignment")

###############################################
#### PLOT 2. Compartment-associated Orders ####
###############################################

# Subset taxa for tree layer
DA_species <- speciesLFC %>% filter(!is.na(compAss)) %$% Species
DA_species.ps <- subset_taxa(moss.ps, Species %in% DA_species)

# table with sample total microbial sequences
relab <- moss.melt %>%
  group_by(Sample) %>% 
  mutate(total = sum(Abundance)) %>% 
  transmute(OTU = OTU,  # keep all these variables
            relAb = (Abundance/total),
            Compart = Compartment,
            Species = Species, 
            Host = Host) %>% 
  inner_join(speciesLFC, by = "Species") %>%  # add compartment association
  filter(!is.na(compAss)) # keep only DA species

# Compute sample-wise mean sequence proportion
relab_comp.df <- relab %>% 
  group_by(Compart, OTU, Species, compAss) %>% 
  summarise(mean_relAb = mean(relAb),
            sd_relAb = sd(relAb)) %>% 
  mutate(across(where(is.character),as.factor)) %>% 
  ungroup %>% arrange(OTU)
  
# PLOT !
compColors <- c('darkgoldenrod4', 'darkolivegreen3')
nOrder <- DA_species.ps@tax_table %>% as.data.frame %$% Order %>% unique %>% length
orderCol <- colorRampPalette(brewer.pal(9, "Set1"))(nOrder+1) # +1 for NAs

tree_taxrank <- function(rank) {
ggtree(DA_species.ps, #layout = "fan", 
       size = 0.2) +
    # Node tips : 
    geom_tippoint(mapping = aes(color = !!sym(rank)), size = 2) +
    scale_colour_manual(values = classCol) +
    #scale_colour_brewer(palette = "Set1") +
    xlim(-1, NA) + # prevent the high-level branches from clustering in the middle
    # Compartment-association tile
    geom_fruit(relab_comp.df,
               geom = geom_tile,
               mapping = aes(y = OTU, x = 0.1, fill = compAss),
               offset = 0, width = 0.1) +
    geom_fruit(relab_comp.df, geom = geom_col, 
                  mapping = aes(y = OTU, x = mean_relAb, fill = Compart),
                  position = position_dodgex(width = 0.7),
                  pwidth = 1, offset = 0.1,
                  # Add a grid behind
                  axis.params=list(
                    axis = "x", text.size = 2, nbreak = 4,
                    text = "Mean relative abundance"
                  ), grid.params=list()) + 
    scale_fill_manual(values = compColors) +
    # geom_fruit(relab_comp.df, geom = geom_errorbar,
    #            mapping = aes(y = OTU, xmin = mean_relAb - sd_relAb,
    #                          xmax =  mean_relAb + sd_relAb),
    #            position = position_dodgex(width = 0.7),
    #            pwidth = 1, offset = 0.1) +
    labs(fill = "Compartment association",
          colour = paste("Bacterial", rank)) +
    theme(plot.margin = margin(t=20, r=20, b=20,l = -200),
          legend.margin = margin(t=30))
}
tree_taxrank('Order')
tree_taxrank('Family')

#########################################################
#### PLOT 3. Acetobacterales distribution across host ####
#########################################################

DA_acetobacterales <- DA_species.ps %>% psmelt %>% 
  filter(Order == "Acetobacterales") %$% Species %>% unique
  
relab %>% 
  filter(Species %in% DA_acetobacterales) %>% 
  ggplot(aes(y = relAb, fill = Host)) +
  geom_boxplot() +
  facet_wrap( ~ Species, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Compartment-associated Acetobacterales relative abundance across host moss.",
       y = "Proportion of sample reads",
       fill = "Host moss species")

##################################################################3

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

# Add a fake taxonomic level to access compAss by species
# tax_table(DA_species.ps) <- 
#   DA_species.ps@tax_table %>% as.data.frame %>% 
#   rownames_to_column("OTU") %>% # otherwise we lose them
#   left_join((relab_comp.df %>% 
#               dplyr::select(Species, compAss) %>% 
#               unique), by = "Species") %>% 
#   column_to_rownames("OTU") %>% as.matrix # tax_table needs a matrix
