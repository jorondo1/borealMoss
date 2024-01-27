library(pacman)
p_load(tidyverse, magrittr, RColorBrewer, phyloseq, colorRamp2, patchwork,
       ggtree, ggtreeExtra)
source("myFunctions.R")
moss.ps <- readRDS("data/R_out/psMossMAGs.RDS")

#####################################
#### PLOT 2. MAGs characteristics ####
#####################################

# DF for MAG quality heatmap 
MAG_names <- moss.ps@tax_table %>% rownames %>% .[grep(".bin.", .)]
hm.mx <- read_tsv("data/novel_quality_scores.txt",
                  col_names = c("MAG", 'comp', 'cont', 'QS')) %>% 
  mutate(MAG = str_remove(MAG, ".fa")) %>% 
  left_join(speciesLFC, by = 'MAG') %>% 
  filter(MAG %in% MAG_names) %>% 
  column_to_rownames("MAG")

# Expand a Brewer palette to 12 colours:
subset.ps <- moss.ps %>% prune_taxa(taxa = MAG_names, .)
nClass <- subset.ps@tax_table %>% as.data.frame %$% Class %>% unique %>% length
classCol <- colorRampPalette(brewer.pal(8, "Set1"))(nClass+1) # +1 for NAs

# Plot the tree
p <- ggtree(subset.ps, layout="fan", 
            size=0.2) +  # Override the colour mapping shape by creating sham geom_point
  xlim(-0.6, NA) + # prevent the high-level branches from clustering in the middle
  geom_tippoint(mapping = aes(color = Class), size = 2.5) +
  scale_colour_manual(values = classCol, na.value = "grey10")

# Add a heatmap
###! !!!!! Check if QS is aligned, there's no anchoring!
gheatmap(p, data = hm.mx["QS"], 
         offset=0.01, width=.1, colnames = FALSE) +
  scale_fill_gradient(low = 'blue4', high = 'gold') +
  #scale_fill_viridis_c(option="E", name="MAG Quality Score") +
  labs(color = "GTDB-Tk Class Assignment")

#########################################
### 3. DIFFERENTIAL ABUNDANCE by HOST ####
#########################################
hostDA <- read_rds('data/R_out/DA_host.RDS')

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
speciesLFC <- readRDS("data/DA_results.RDS") %>% 
  transmute(LFC = lfc_CompartmentGreen, 
            Species = taxon, 
            compAss = case_when(lfc_CompartmentGreen>0 ~ "Green",
                                lfc_CompartmentGreen<0 ~ "Brown")) %>% 
  right_join(moss.ps %>% # identifier \ species association table
               tax_table %>% as.data.frame %>% 
               select(Species) %>% rownames_to_column("MAG"),
             by = 'Species')

# Subset taxa for tree layer
DA_species <- speciesLFC %>% filter(!is.na(compAss)) %$% Species
DA_species.ps <- subset_taxa(moss.ps, 
                             Species %in% DA_species)

rank <- "Order"
n <- DA_species.ps@tax_table %>% as.data.frame %>% .[rank] %>% unique %>% dim %>% .[1]

### Taxonomic tree (generated first to establish species factor levels)
tree.p <- ggtree(DA_species.ps,size = 0.2) +
  geom_tippoint(mapping = aes(color = !!sym(rank)), size = 3) +
  #geom_tiplab(align=TRUE) +
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
  filter(!is.na(label)) %$% Species# NAs at non-integer positions?!

### Waterfall plot :

DA_comp.df <- read_rds("data/R_out/DA_comp.RDS") %>% 
  mutate(Species = factor(taxon, levels = orderedSpecies), .keep = 'unused')
  
# Subset the ~80% largest effect:
# subset <- DA_comp.df %>% 
#   transmute(Species = Species, 
#             lfc = abs(lfc_CompartmentGreen)) %>%
#   arrange(desc(lfc)) %>% 
#   head(round((DA_comp.df %>% dim %>% .[1])*0.7)) %$% Species

sfx <- "CompartmentGreen"
waterfall.p <- DA_comp.df %>% 
  #filter(Species %in% subset) %>% 
  ggplot(aes(y = Species, x = !!sym(paste0("lfc_",sfx)), fill = direct)) + 
  geom_bar(stat = "identity", width = 1, color = "white",
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(xmin = !!sym(paste0("lfc_",sfx)) - !!sym(paste0("se_",sfx)), 
                    xmax = !!sym(paste0("lfc_",sfx)) + !!sym(paste0("se_",sfx))), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(y = NULL, x = "Log fold change") + 
  scale_fill_manual(values = c("Positive LFC" = compColours[2], "Negative LFC" = compColours[1]))+
  scale_color_discrete(name = NULL) +
  #scale_y_discrete(position = 'right') +
  theme_minimal() + 
  theme(panel.grid.minor.y = element_blank(),
        plot.margin = margin(t = 10, r = 10, b = 10, l = -5, unit = "pt")) +
  guides(fill = FALSE)

tree.p + waterfall.p +
  plot_layout(#guides = "collect",
              design = "AAAABBB") +
  plot_annotation(caption = 'Significantly abundant at p<0.01 (ajdusted).\nRestricted to species with >10% relative abundance that passed the sensitivity analysis.')
