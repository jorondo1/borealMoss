library(pacman)
p_load(ape, tidyverse, magrittr, ggtree, phyloseq,
       RColorBrewer, ggtreeExtra, ggnewscale)
source("myFunctions.R")


###############################################
#### PLOT 2. Compartment-associated Orders ####
###############################################
moss.melt <- readRDS("data/psMossMAGs.RDS") %>% psmelt

# Subset taxa for tree layer
DA_species <- speciesLFC %>% filter(!is.na(compAss)) %$% Species
DA_species.ps <- subset_taxa(moss.ps, Species %in% DA_species)

# table with sample total microbial sequences
relab <- moss.melt %>%
  group_by(Sample) %>% 
  mutate(total = sum(Abundance)) %>% 
  mutate(relAb = (Abundance/total),
            Compart = Compartment, .keep = 'unused') %>% 
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
tree_taxrank <- function(ps, rank) {
  n <- ps@tax_table %>% as.data.frame %>% .[rank] %>% unique %>% dim %>% .[1]

  ggtree(DA_species.ps,size = 0.2) +
    # Node tips : 
    geom_tippoint(mapping = aes(color = !!sym(rank)), size = 2) +
    geom_tiplab(as_ylab = TRUE) + # add labels to the right
    scale_colour_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(n) ) +
    xlim(-1, NA) + # horizontal tree aligment
    # Compartment-association tile
    geom_fruit(relab_comp.df,
               geom = geom_tile,
               mapping = aes(y = OTU, x = 0.1, fill = compAss),
               offset = 0, width = 0.1) +
    # Proportion of reads barplot
    geom_fruit(relab_comp.df, geom = geom_col, 
               mapping = aes(y = OTU, x = mean_relAb, fill = Compart),
               position = position_dodgex(width = 0.7),
               pwidth = 1, offset = 0.1,
               # add a grid
               grid.params=list()) + 
    #geom_fruit() +
    scale_fill_manual(values = compColours) +
    labs(fill = "Compartment",
         colour = paste("Bacterial", rank),
         x = 'Mean proportion of sample sequences') +
    theme(plot.margin = margin(t=20, r=20, b=20,l = -200),
          legend.margin = margin(t=30),
          axis.title.x = element_text(hjust = 0.95))
}
tree_taxrank(DA_species.ps, 'Order')
tree_taxrank(DA_species.ps, 'Family')
tree_taxrank(DA_species.ps, 'Class')

#########################################################
#### PLOT 3. Acetobacterales distribution across host ####
#########################################################

clade_host <- function(rank, taxa) {
  DA <- DA_species.ps %>% psmelt %>% 
    filter(!!sym(rank) == taxa) %$% Species %>% unique
  
  relab %>% 
    filter(Species %in% DA &
             Compart == compAss) %>%
    ggplot(aes(y = relAb, fill = Host)) +
    geom_boxplot() +
    #facet_grid(rows = 'Species') +
    facet_wrap( ~ Species, scales = "free") +
    theme_minimal() +
    theme(axis.text.x = element_blank()) +
    labs(title = paste("Compartment-associated", taxa, "relative abundance across host moss."),
         y = "Proportion of sample reads",
         fill = "Host moss species")
}

clade_host("Order", "Acetobacterales")
clade_host("Family", "Beijerinckiaceae")
clade_host("Family", "Xanthobacteraceae")

#########################################################
#### PLOT 3. Compartment association across host for selected families ####
#########################################################
families <- c("Beijerinckiaceae", "Xanthobacteraceae")

clade_host_compart <- function(clade, rank, taxa) {
  relab %>% 
    filter(!!sym(clade) %in% taxa) %>% 
    group_by(Sample, !!sym(rank), Host, Compart) %>% 
    summarise(relAb = sum(relAb)) %>% 
    # PLOT 
    ggplot(aes(x = Host, y = relAb)) + 
    geom_boxplot(aes(fill = Compart)) + 
    geom_point(aes(colour=Compart), size = 0.7,
               position = position_jitterdodge()) +
    scale_fill_manual(values = compColours) +
    scale_colour_manual(values = c('black', 'black')) +
    facet_grid(rows = paste0(rank,' ~ Host'), scales = "free") +
    theme_light() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    labs(y = "Proportion of sample sequences counts")
}
clade_host_compart('Family', 'Family', families) +
  labs(title = "Clade sequence count distribution by Host and Family")
clade_host_compart('Family', 'Genus', 'Acetobacteraceae') +
  labs(title = "Acetobacteraceae sequence count distribution by Host and Genus")
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
