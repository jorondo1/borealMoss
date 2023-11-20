library(pacman)
p_load(tidyverse, magrittr, DESeq2, vegan, RColorBrewer, bestNormalize)
source("myFunctions.R")

psMossGTDB <- readRDS("data/psMossGTDB.RDS")
psMossMAGs <- readRDS("data/psMossMAGs.RDS")

MAGs_melt <- psMossMAGs %>% psmelt %>% 
  select(OTU, Sample, Abundance, Compartment, Microsite, Host, Domain:Species)
######################
### ALPHA DIVERSITY ###
######################

# Shannon diversity calculated by rarefaction
# This function normalises if the Shannon diversity is not normal

rare_MAGs <- rarefy_even_depth(psMossMAGs,
                                    #  sample.size = 260000, 
                                    replace = FALSE,
                                    verbose = TRUE,  
                                    trimOTUs = TRUE,
                                    rngseed = 4466) %>% 
  estimate_richness(measure = c("Shannon")) 
shapiro.test(rare_MAGs$Shannon)

shannon.norm <- bestNormalize(rare_MAGs$Shannon) #orderNorm
shapiro.test(shannon.norm$x.t)

div_plots <- data.frame(sample_data(psMossMAGs),
           Shannon=shannon.norm$x.t)

cbind(sample_data(psMossMAGs), shannon.norm$x.t) %>% 
  mutate(Compartment = as.factor(Compartment)) 

div_plots %>% 
  ggplot(aes(x = Compartment, y = Shannon)) +
  geom_boxplot(aes(fill = Compartment)) + 
  facet_wrap("Host", ncol = 4, 
             labeller = label_wrap_gen(width=10)) +
  theme_bw() +
  ggtitle("TAXONOMIC DIVERSITY OF MOSS COMPARTMENTS") + 
  ylab("SHANNON ALPHA-DIVERSITY") + 
  theme(axis.text.x  = element_text (size=12, color="black", angle=45, vjust = 1, hjust = 1),
        axis.text.y  = element_text (size=12, color="black"),
        axis.title.x  = element_text(size=12, color="black"),
        axis.title.y  = element_text(size=12, color="black"),
        strip.text = element_text(size = 14)) + 
  scale_fill_manual(values = c("tan4", "springgreen4")) + 
  guides(fill="none") + coord_fixed(2)

#####################
### BETA DIVERSITY ###
#####################

vst.mx <- psMossMAGs %>% 
  phyloseq_to_deseq2(~Compartment) %>% # DESeq2 object
  estimateSizeFactors(., geoMeans = apply(
    counts(.), 1, function(x) exp(sum(log(x[x>0]))/length(x)))) %>% 
  DESeq2::varianceStabilizingTransformation(blind=T) %>% # VST
  SummarizedExperiment::assay(.) %>% t %>% 
  { .[. < 0] <- 0; . }

dist.mx <- vegan::vegdist(vst.mx, distance = "jaccard")
PCoA <- capscale(dist.mx~1, distance = "jaccard")

# Plot data
plot.df <- data.frame(PCOA1 = PCoA %>% scores %$% sites %>% .[,1], 
                      PCOA2 = PCoA %>% scores %$% sites %>% .[,2]) %>%
  cbind(subset.ps %>% sample_data %>% data.frame)

PCoA$CA$eig[1:3]/sum(PCoA$CA$eig)
eig.test <- PCoA$CA$eig
eig.test[1]/sum(abs(eig.test)) 
eig.test[2]/sum(abs(eig.test))
eig.test[3]/sum(abs(eig.test))

div_plots$MDS1<-scores(PCoA)$sites[,1]
div_plots$MDS2<-scores(PCoA)$sites[,2]

# Ordination plot between compartments
div_plots %>% 
  ggplot(aes(MDS1, MDS2, colour = Compartment)) + 
  stat_ellipse(level=0.9, geom = "polygon", alpha = 0.18, aes(fill = Compartment)) +   
  geom_point(size = 5, aes(shape = Host)) + 
  ggtitle("TAXONOMIC BETA-DIVERSITY BETWEEN SITE TYPES") +
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 14)) + 
  guides(fill="none") + 
  scale_fill_manual(values = c("tan4", "springgreen3")) + 
  scale_color_manual(values = c("lightsalmon4", "springgreen4")) +  
  coord_fixed()

# Ordination plot between HOST SPECIES
div_plots %>% 
  ggplot(aes(MDS1, MDS2, colour = Host)) + 
  stat_ellipse(level=0.9, geom = "polygon", alpha = 0.18, aes(fill = Host)) +   
  geom_point(size = 5, aes(shape = Compartment)) + 
  ggtitle("TAXONOMIC BETA-DIVERSITY BETWEEN SITE TYPES") +
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 14)) + 
  guides(fill="none") + 
#  scale_fill_manual(values = c("tan4", "springgreen3")) + 
#  scale_color_manual(values = c("lightsalmon4", "springgreen4")) +  
  coord_fixed()

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

RCLR <- trans.fun(psMossMAGs,'rclr',topSpecies)

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
