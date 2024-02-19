library(pacman)
p_load(ape, tidyverse, magrittr, RColorBrewer, colorRamp2, patchwork,
       ggtree, ggtreeExtra, treeio, ggnewscale, cowplot)
source("myFunctions.R")
moss.ps <- readRDS("data/R_out/mossMAGs.RDS")

###################################
#### PLOT 1. Community Overview ####
###################################

topTaxa <- function(psmelt, comp, taxLvl, topN) {
# Compute top taxa, aggregate residual in others
  psmelt %>% 
    filter(Compartment == comp) %>% # subset compartment
    group_by(!!sym(taxLvl)) %>% # group by tax level
    summarise(relAb = mean(relAb)) %>% # find top abundant taxa
    arrange(desc(relAb)) %>% 
    mutate(aggTaxo = as.factor(case_when( # aggTaxo will become the plot legend
      row_number() <= topN ~ !!sym(taxLvl), #+++ We'll need to manually order the species!
      row_number() > topN ~ 'Others'))) # +1 to include the Others section!
}

df_comm <- function(psmelt, comp, taxLvl, topTaxa) {
    left_join(psmelt, topTaxa, by=taxLvl) %>% # use topTaxa to join aggTaxo variable
    filter(Compartment == comp) %>% 
    aggregate(Abundance~aggTaxo+Host+Compartment, # Abundance is aggregated...
              data=., FUN = sum) %>% # ...sums "Others" taxa by specified variables
    mutate(aggTaxo = factor(aggTaxo,# reorder 
                            levels = topTaxaLvls)) 
}

# Desired taxonomic aggregation level :
taxLvl <- 'Class'

# Preliminary dataset with variables of interest
MAGs_melt <- moss.ps %>% psmelt %>%
  select(Sample, Abundance, Compartment, Host, Domain:Species) %>% 
  # Compute Class level abundance :
  group_by(Sample, across(all_of(taxLvl))) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), 
            .groups = 'drop') %>% 
  group_by(Sample) %>% # Convert to relative abundance
  mutate(relAb = Abundance/sum(Abundance)) %>% ungroup

topN=9
topTaxa_Brown <- topTaxa(MAGs_melt, 'Brown', taxLvl, topN)
topTaxa_Green <- topTaxa(MAGs_melt, 'Green', taxLvl, topN)

# Create ordered list of taxa
topTaxaLvls <- rbind(topTaxa_Brown, topTaxa_Green) %>% 
  group_by(aggTaxo) %>% 
  aggregate(relAb ~ aggTaxo, data = ., FUN = sum) %>% 
  arrange(relAb) %$% aggTaxo %>% as.character 

df_Brown <- df_comm(MAGs_melt, 'Brown', taxLvl, topTaxa_Brown)
df_Green <- df_comm(MAGs_melt, 'Green', taxLvl, topTaxa_Green)

# Plot !
mycolors1 <- colorRampPalette(brewer.pal(6, "Set3"))(topN+1)
ggplot(df_Green) +
  labs(title = 'Green section') +
ggplot(df_Brown) +
  labs(title = 'Brown section') +
patchwork::plot_layout(guides = 'collect', 
                       axes = 'collect',
                       ncol = 1) &
  geom_bar(stat = "identity", 
           position = "fill", 
           aes(x = Host, y = Abundance, fill = aggTaxo)) &
  labs(fill = taxLvl, 
       y = "Proportion of k-mer counts", 
       x = 'Host moss species')  &
  scale_fill_manual(values = mycolors1, breaks = topTaxaLvls) &
  scale_x_discrete(labels = (c(expression(italic("D. undulatum")),
                               expression(italic("P. commune")),
                               expression(italic("P. juniperinum")),
                               expression(italic("P. piliferum"))))) &
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

#####################################
#### PLOT 2. MAGs characteristics ####
#####################################

# extract MAG names
MAG_names <- moss.ps@tax_table %>% rownames %>% .[grep(".bin.", .)] 

# Extract taxonomy
taxLabels <- moss.ps@tax_table %>% as.data.frame %>% 
  rownames_to_column("label") %>% 
  select(label, Domain:Species) %>% 
  tibble

MAG.tree <- read.tree("data/RAxML_bestTree.MAGs_refined.tre")
# Subset tree to MAGs only and add taxonomy
sub.tree <- MAG.tree %>% 
  drop.tip(setdiff(MAG.tree$tip.label, MAG_names)) %>% 
  as_tibble %>% 
  full_join(taxLabels %>% filter(label %in% MAG_names), by = 'label') %>% 
  as.treedata

# Expand a Brewer palette to 12 colours:
nClass <- moss.ps %>% prune_taxa(taxa = MAG_names, .) %>% 
  tax_table %>% as.data.frame %$% Class %>% unique %>% length
classCol <- colorRampPalette(brewer.pal(9, "Set1"))(nClass)

# Plot the tree
p <- ggtree(sub.tree, layout="fan", size=0.1) +
  xlim(-0.2, NA) + # prevent the high-level branches from clustering in the middle
  geom_tippoint(mapping = aes(color = Class), size = 2.5) +
  scale_colour_manual(values = classCol) 
# add MAG data:
hm.mx <- read_tsv("data/R_out/MAG_summary.tsv") %>% 
  # only keep MAGs that are in our original dataset
  filter(MAG %in% MAG_names) 

# Add a heatmap
p + geom_fruit(data = hm.mx, geom = geom_tile,
                 mapping = aes(y = MAG, fill = MBP),
                 offset = 0.1, width = 0.1) +
  scale_fill_gradient(low = "pink1", high = 'purple4', name="Length (Mbp)") +
  
  new_scale_fill() + 
  geom_fruit(data = hm.mx, geom = geom_tile,
               mapping=aes(y=MAG, fill = GC),
               offset = 0.1,width = 0.1) +
  scale_fill_gradient(low = 'darkgreen', high = 'gold', name="%GC Content") +
  
  new_scale_fill() +
  geom_fruit(data = hm.mx, geom = geom_tile,
             mapping = aes(y = MAG, fill = `L50 (kb)`),
             offset = 0.1, width = 0.1) +
  scale_fill_gradient(low = 'lightblue', high = 'darkblue') +
  
  new_scale_fill() +
  geom_fruit(data = hm.mx, geom = geom_tile,
             mapping = aes(y = MAG, fill = QS),
             offset = 0.1, width = 0.1) + 
  scale_fill_gradient(low = "darkred", high = "salmon1", name="Quality Score")
  
# theme(legend.position = 'bottom',
#       legend.direction = 'vertical')

#########################################
### 3. DIFFERENTIAL ABUNDANCE by HOST ####
#########################################
hostDA <- read_rds('data/R_out/DA_host_results.RDS')

# Which tax level are we showing in the legend?
taxLvl <- "Class"

# sort Species by taxLvl (descending so species are top-to-bottom on y axis)
speciesLvl <- hostDA %>% arrange(desc(!!sym(taxLvl)), taxon) %$% taxon %>% unique

DA_host.df <- hostDA %>% 
  #filter(Group %in% c("P_commune", "P_juniperinum", "P_piliferum")) %>% 
  # remove all orders for which diff is FALSE in every group :
  group_by(taxon) %>% 
  filter(!all(diff == FALSE)) %>% ungroup() %>% 
  # reorder taxa by taxLvl
  mutate(taxon = factor(taxon, levels = speciesLvl))

# MAIN DAA PLOT :
DA.p1 <- DA_host.df %>% 
  ggplot(aes(x = Group, y = taxon, fill = lfc)) +
  geom_tile() +
  scale_fill_gradient2(low = "#ff7f0e", high = "#1f77b4", mid = "grey95", 
                       midpoint = 0) +
  geom_text(aes(Group, taxon, label = round(lfc,2), color=textcolour)) +
  scale_color_identity(guide = FALSE) + 
  theme_void() + 
  theme(axis.text.x = element_text(margin = margin(t = 5, r = 5, b = 5, l = 5))) +
  labs(#title = paste0("Differential abundance analysis at the ",taxRank," level."),
    x = "", y = "",
    fill = "Log-fold change\nin abundance\nrelative to\nD. dicranum.")

# taxLvl -coloured tile:
DA.p2 <- DA_host.df %>% 
  ggplot(aes(x = '1', y = taxon, fill = !!sym(taxLvl))) +
  geom_tile() + theme_void() + 
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y = element_text(hjust = 1)) +
  scale_fill_brewer(palette = "Paired")

DA.p2 + DA.p1 + 
  plot_layout(
    guides = "collect",
    design = "ABBBBBBBBBBBB")

################################################
### 4. DIFFERENTIAL ABUNDANCE by COMPARTMENT ####
################################################

# Add LFC (from DAA) to significant species
speciesLFC <- readRDS("data/R_out/speciesLFC_comp.RDS")
tree <- read.tree("data/RAxML_bestTree.genomes_refined.tre") 
# Subset taxa for tree layer
DA_species <- speciesLFC %$% MAG
DA_sub.tree <- tree %>% 
  drop.tip(setdiff(tree$tip.label, DA_species)) %>% 
  as_tibble %>% 
  # add taxonomy :
  full_join(taxLabels %>% filter(label %in% DA_species), by = 'label') %>% 
  as.treedata # because.

rank <- "Class"
n <- DA_sub.tree@data %>% as.data.frame %>% .[rank] %>% unique %>% dim %>% .[1]

### Taxonomic tree (generated first to establish species factor levels)
tree.p <- ggtree(DA_sub.tree,size = 0.2) 

orderedRank <- tree.p$data %>% filter(!is.na(label)) %>% 
  arrange(desc(y)) %>% select(!!sym(rank)) %>% 
  # dplyr way of collapsing a tibble into a vector is pull() :
  pull %>% unique

tree.p$data %<>% mutate(!!sym(rank) := # dynamic management of variable name
                         factor(!!sym(rank), 
                                levels = orderedRank))

tree.p <- tree.p +
  geom_tippoint(mapping = aes(color = !!sym(rank)), size = 3) +
  #geom_tiplab(align=TRUE, aes(label = Species)) +
  scale_colour_manual(values = colorRampPalette(brewer.pal(9, "Set1"))(n) ) +
  scale_fill_manual(values = compColours) +
  labs(fill = "Compartment",
       colour = paste("Bacterial", rank)) +
  theme(plot.margin = margin(t=20, r=0, b=20,l = 20),
        legend.margin = margin(t=30),
        axis.title.x = element_text(hjust = 0.95)) +
    theme(legend.position = c(.1,.8))

# Extract species levels for alignment of next plots :
orderedSpecies <- tree.p$data %>% select(y, label, Species) %>% 
  unique %>% arrange(y) %>% # sort by plot position
  filter(!is.na(label)) %$% Species # NAs at non-integer positions?!

### Waterfall plot :
speciesLFC %<>% mutate(Species = factor(Species, levels = orderedSpecies))
  
waterfall.p <- speciesLFC %>% 
  #filter(Species %in% subset) %>% 
  ggplot(aes(y = Species, x = LFC, fill = compAss)) + 
  geom_bar(stat = "identity", width = 1, color = "white",
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(xmin = LFC - SD, 
                    xmax = LFC + SD), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(y = NULL, x = "Log fold change", fill = 'Moss section\nassociation') + 
  scale_fill_manual(values = c("Green" = compColours[2], "Brown" = compColours[1]))+
  scale_color_discrete(name = NULL) +
  theme_minimal() + 
  theme(panel.grid.minor.y = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = -5, unit = "pt"),
        legend.position = c(.8,.95),
        legend.background = element_rect(colour='black', fill='white', linewidth=0.2),
        legend.title = element_text(size = rel(.8)))

tree.p + waterfall.p +
  plot_layout(#guides = "collect",
              design = "AAAABBB") +
  plot_annotation(caption = 'Significantly abundant at p<0.05 (ajdusted).\nRestricted to species with >10% relative abundance that passed the sensitivity analysis.')

