library(pacman)
p_load(tidyverse, magrittr, DESeq2, vegan, RColorBrewer, microbiome,
       bestNormalize, patchwork)
source("myFunctions.R")

mossGTDB.ps <- readRDS("data/R_out/mossGTDB.RDS")
mossMAGs.ps <- readRDS("data/R_out/mossMAGs.RDS")

MAGs_melt <- mossMAGs.ps %>% psmelt %>% 
  select(OTU, Sample, Abundance, Compartment, Microsite, Host, Domain:Species)

######################
### ALPHA DIVERSITY ###
######################

# Shannon diversity calculated by rarefaction
# This function normalises if the Shannon diversity is not normal
div.boxplot <- function(ps, title) {
  rare_MAGs <- rarefy_even_depth(ps,
                                      #  sample.size = 260000, 
                                      replace = FALSE,
                                      verbose = TRUE,  
                                      trimOTUs = TRUE,
                                      rngseed = 4466) %>% 
    estimate_richness(measure = c("Shannon")) 
  shapiro.test(rare_MAGs$Shannon)
  
  shannon.norm <- bestNormalize(rare_MAGs$Shannon) #orderNorm
  shapiro.test(shannon.norm$x.t)
  
  div_plots <- data.frame(sample_data(ps),
             Shannon=shannon.norm$x.t)
  
  cbind(sample_data(ps), shannon.norm$x.t) %>% 
    mutate(Compartment = as.factor(Compartment))
  
  div_plots %>% 
    ggplot(aes(x = Compartment, y = Shannon)) +
      geom_boxplot(aes(fill = Compartment)) +
      facet_wrap("Host", ncol = 4) +
      theme_bw() +
      ggtitle(title) +
      ylab("Shannon Alpha-diversity") +
      # theme(axis.text.x  = element_text (size=12, color="black", angle=45, vjust = 1, hjust = 1),
      #       axis.text.y  = element_text (size=12, color="black"),
      #       axis.title.x  = element_text(size=12, color="black"),
      #       axis.title.y  = element_text(size=12, color="black"),
      #       strip.text = element_text(size = 14)) +
      scale_fill_manual(values = c("tan4", "springgreen4")) +
      guides(fill="none") + coord_fixed(2)
}

alpha_MAGs.plot <- div.boxplot(mossMAGs.ps, "Alpha diversity with nMAGs.")
alpha_GTDB.plot <- div.boxplot(mossGTDB.ps, "Alpha diversity without nMAGs.")

alpha_MAGs.plot + alpha_GTDB.plot
  
#####################
### BETA DIVERSITY ###
#####################

plotOrd.fun <- function(ps, title, var, dist) {
    # PCoA custom function :
  df <- pcoa.fun(ps, var, vst.fun(ps, 'Compartment'), dist) 
  
  # Ordination plot between compartments :
  df$metadata %>% 
    ggplot(aes(x = PCoA1, y = PCoA2, colour = !!sym(var))) + 
    stat_ellipse(level=0.9, geom = "polygon", alpha = 0.18, aes(fill = !!sym(var))) +   
    geom_point(size = 5, aes(shape = Host)) + 
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(size = 18),
          legend.title = element_text(colour="black", size=16, face="bold"),
          legend.text = element_text(colour="black", size = 14)) + 
    guides(fill="none") + 
    labs(
      x = paste0(
      "PCoA 1 [", round(100*df[['eig']][1]/sum(abs(df[['eig']])),1), "% ]"),
      y = paste0(
        "PCoA 2 [",round(100*df[['eig']][2]/sum(abs(df[['eig']])),1),"% ]")) +
    scale_fill_manual(values = compColours) + 
    scale_color_manual(values = compColours) 
}
# Variance-stabilizing transformation custom function :

beta_MAGs.plot <- plotOrd.fun(mossMAGs.ps, "Beta diversity with nMAGs.", 'Compartment', 'bray')
beta_GTDB.plot <- plotOrd.fun(mossGTDB.ps, "Beta diversity without nMAGs.", 'Compartment', 'bray')

beta_MAGs.plot + beta_GTDB.plot + plot_layout(guides = 'collect')

### Compare with or without D. undulatum:
beta_MAGs_all.plot <- ord.fun(mossMAGs.ps, "Beta diversity with Dicranum undulatum")
mossMAGs_poly.ps <-  mossMAGs.ps %>% 
  # remove D_undulatum samples:
  prune_samples(sample_data(.)$Host != "D_undulatum",.) %>%
  # remove Taxa absent from all samples 
  prune_taxa(taxa_sums(.) > 0,.) 

beta_MAGs_Polytrichum.plot <- mossMAGs_poly.ps %>% 
#  prune_samples(sample_data(.)$Host != "P_commune",.) %>% 
#  prune_samples(sample_data(.)$Location != "Nemaska community intersection",.) %>% 
#  prune_samples(sample_data(.)$Location != "Chemin PK",.) %>% 
  # ordinate & plot
  ord.fun("Beta diversity without Dicranum undulatum")

beta_MAGs_all.plot + beta_MAGs_Polytrichum.plot +
  guides(shape = FALSE) +
  plot_layout(guides = 'collect')

# Ordination plot between HOST SPECIES
df$metadata %>% ##<<<<!!! this df is redefined every time ord.fun is executed!
  ggplot(aes(PCoA1, PCoA2, colour = Host)) + 
  stat_ellipse(level=0.9, geom = "polygon", alpha = 0.18, aes(fill = Host)) +   
  geom_point(size = 5, aes(shape = Compartment)) + 
  ggtitle("Beta-diversity across host moss species") +
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 14)) + 
  guides(fill="none") + 
  coord_fixed()

# Ordination vs temp gradients
df <- vst.fun(mossMAGs_poly.ps, 'Compartment') %>% 
  pcoa.fun(mossMAGs_poly.ps, 'Compartment', ., 'bray')

df$metadata %>% 
  ggplot(aes(PCoA1, PCoA2, colour = Canopy)) + 
  geom_point(size = 5, aes(shape = Compartment)) + 
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 14)) + 
  guides(fill="none") + 
  coord_fixed()

# PERMANOVA
h1 <- with(df$metadata, how(nperm=9999, blocks = Location)) 
adonis2(df$dissMx ~ SoilpH + SoilTemp + Canopy + Compartment*Host, 
        data = df$metadata,
        permutations=h1, method="bray")

# Aitchison way
aitchison.df <- mossMAGs_poly.ps %>% 
  otu_table %>% t %>% # transpose otherwise it computes distance between species
  pcoa.fun(mossMAGs_poly.ps, 'Compartment', ., 'robust.aitchison')

aitchison.df$metadata %>% 
  ggplot(aes(PCoA1, PCoA2, colour = Canopy)) + 
  geom_point(size = 5, aes(shape = Compartment)) + 
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 14)) + 
  guides(fill="none") + 
  coord_fixed()

h1 <- with(aitchison.df$metadata, how(nperm=9999, blocks = Location)) 
adonis2(aitchison.df$dissMx ~ SoilpH + SoilTemp + Canopy + Compartment*Host, 
        data = aitchison.df$metadata,
        permutations=h1, method="robust.aitchison")


##########################
### COMMUNITY OVERVIEW ###
##########################

# Just a check!
MAGs_melt %>% 
  mutate(Group = ifelse(endsWith(Species, "MAG"), "MAG", "Default")) %>%
  group_by(Group) %>%
  summarise(TotalAbundance = sum(Abundance))
# There are about 4x more sequences attributed to our novel MAGs than there 
# are to all 85k species representatives of the GTDB taxonomy

### We find the top taxa and assign the rest to "Others"
comm_barplot <- function(melt, comp, taxLvl, topN) {
    
  mycolors1 <- colorRampPalette(brewer.pal(8, "Set2"))(topN+1)
  
  topTaxa <- melt %>%  
    group_by(!!sym(taxLvl)) %>% # melt the table and group by tax level
    summarise(Abundance = sum(Abundance)) %>% # find most overall abundant taxa
    arrange(desc(Abundance)) %>%  # Species them 
    mutate(aggTaxo = as.factor(case_when( # aggTaxo will become the plot legend
      row_number()<=topN ~ !!sym(taxLvl), #+++ We'll need to manually order the species!
      row_number()>topN~'Others'))) %>%  # +1 to include the Others section!
    select(-Abundance)
  
  topTaxaLvls <- topTaxa$aggTaxo %>% head(topN) %>% as.character %>% sort %>% c("Others")
  
  ### We create an object that will be fed into ggplot
  df_comm <- melt %>% 
    left_join(.,topTaxa, by=taxLvl) %>% # use topTaxa to join aggTaxo variable
    aggregate(Abundance~aggTaxo+Host+Compartment, # Abundance is aggregated...
              data=., FUN = sum) %>% # ...sums "Others" taxa by specified variables
    mutate(aggTaxo = factor(aggTaxo,# reorder 
                              levels = topTaxaLvls))
  
  df_comm %>% ggplot(aes(x = Host, y = Abundance, fill = aggTaxo)) +
    geom_bar(stat="identity", position="fill") + # change to "stack" to plot total counts
    labs(title=paste0('Moss microbiome microbial ', taxLvl, ' in ',comp, ' compartment.'), 
         fill=taxLvl, 
         x="Moss Species", 
         y="Proportion of k-mer counts") +
    # facet_grid(rows = vars(Compartment),# RHS of the tilde is columns
    #            scales="free", # hide samples without counts
    #            space="free_x",
    #            switch="both") + # constant bar width despite different # samples per group
    # some formatting...
    scale_fill_manual(values = mycolors1,
                      breaks = topTaxaLvls) +
    theme_light(base_family="Baskerville", base_size = 20)+
    theme(axis.line = element_blank(), 
          axis.ticks = element_blank(),
          plot.margin = unit(c(1,0.5,1,2), "cm"), 
          plot.caption.position="plot",
          axis.title.y.left = element_text(vjust = 3), 
          legend.spacing.y = unit(0.7, 'cm')) +
    guides(fill = guide_legend(byrow = TRUE))
}

comm_barplot(MAGs_melt, 'Brown', 'Family', 14)
comm_barplot(MAGs_melt, 'Green', 'Family', 14)
comm_barplot(MAGs_melt, 'Brown', 'Order', 10)
comm_barplot(MAGs_melt, 'Green', 'Order', 10)

# Which MAGs are cyanobacteria? 
MAGs_melt %>% filter(Phylum == 'Cyanobacteria') %>% 
  select(OTU, Family, Genus, Species, Compartment) %>% unique



# Top species per host and compartment 
MAGs_melt %>% 
  dplyr::select(Host, Abundance, Compartment, Species, OTU) %>% 
  group_by(Host, Compartment, Species) %>% 
  summarise(Abundance = sum(Abundance)) %>% ungroup %>% 
  group_by(Host, Compartment) %>% 
  mutate(TotalAbundance = sum(Abundance)) %>% 
  arrange(desc(Abundance)) %>%
  mutate(Rank = row_number()) %>% 
  filter(Rank <= 3) %>% 
  mutate(PercentageOfTotal = (Abundance / TotalAbundance) * 100) %>% 
  select(-Abundance, -TotalAbundance) %>% 
  arrange(Host, Compartment) %>% View

###########################
#### HEATMAP ##############
###########################

# Top secies per sample ###
topSpecies <- MAGs_melt %>% 
  group_by(Sample) %>%
  arrange(desc(Abundance)) %>%
  slice_head(n = 1) %>% 
  ungroup() %$% OTU %>% unique # Use OTU because two OTUs may have the 
            # same Species field if they are unknown species of the same genus

RCLR <- trans.fun(mossMAGs.ps,'rclr',topSpecies)

heatmap_list <- list(); for (host in unique(RCLR@sam_data$Host)) {
  for (comp in unique(RCLR@sam_data$Compartment)) {
    
    subset_data <- subset_samples(RCLR, Compartment == comp & Host == host)
    
    # Create a name for the heatmap combining host and compartment info
    heatmap_name <- paste("Host:", host, "- Compartment:", comp)
    
    hm <- subset_data %>% 
      comp_heatmap(show_row_dend = F,
                   show_column_dend = F,
                   name = heatmap_name,
                   show_heatmap_legend = FALSE) # Add the name parameter here
    
    # Store each heatmap in the list
    heatmap_list[[paste0(host, comp, "_HMap")]] <- hm
  }
}

# Combine all heatmaps into a single plot
if (length(heatmap_list) > 0) {
  # Concatenate heatmaps using the + operator
  combined_heatmap <- Reduce(`+`, heatmap_list)
  # Draw the combined heatmap
  draw(combined_heatmap, heatmap_legend_side = "bot")
}
