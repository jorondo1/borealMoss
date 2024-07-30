### AUTHOR : SARAH ISHAK ######################################################

library(pacman)
p_load(tidyverse, phyloseq, DESeq2, vegan, RColorBrewer, bestNormalize, 
       wesanderson, sjPlot, car, lme4, emmeans, DHARMa, reshape2, 
       rstatix, ggpubr, circlize, picante, ggeffects,ComplexHeatmap)
source("scripts/myFunctions.R")

#### Loading data ####
psMossMAGs <- readRDS("data/R_out/mossMAGs.RDS")
asv <- as.data.frame(psMossMAGs@otu_table)

# Remove species with fewer than 100 counts total 
length(which(apply(asv,1,sum)<100))
asv100 <- asv[apply(asv,1,sum)>=100,]  # sum of sequences per ASVs
hist(log10(apply(asv,1,sum))) 
hist(log10(apply(asv100,1,sum))) 
hist(log10(apply(asv100,2,sum)))

# New phyloseq 
metadata <- as.data.frame(psMossMAGs@sam_data)
taxa <- as.data.frame(psMossMAGs@tax_table)

ps100 <- phyloseq(sample_data(metadata),
                     otu_table(asv100, taxa_are_rows = TRUE), 
                     tax_table(as.matrix(taxa)))

#########################
#### ALPHA DIVERSITY #### 
########################

# Shannon diversity calculated by rarefaction
rare_MAGs <- rarefy_even_depth(ps100,sample.size = 24900, replace = FALSE, verbose = TRUE,  trimOTUs = TRUE, rngseed = 4466)
rare_rich <- estimate_richness(rare_MAGs, measure = c("Shannon")) 

samsums_rare <- sample_sums(rare_MAGs)
shapiro.test(rare_rich$Shannon)
# the shapiro test came out normal for Shannon (0.024)...so we wouldn't need to normalise it...

palette <- wes_palette("AsteroidCity1", n=4, type = "discrete")  
# modelling ####
div_plots <- data.frame(sample_data(rare_MAGs),
                        Shannon=rare_rich$Shannon) %>% 
  mutate(Moss = Host) %>% 
  separate(Host, c("Genus", "Species"), "_") %>% 
  mutate(
    Genus = case_when(
      Genus == "P" ~ "Polytrichum",
      Genus == "D" ~ "Dicranum",
      TRUE ~ Genus
    )
  ) %>% dplyr::rename(Host = Moss) %>% 
  mutate(across(where(is.character), factor))

# full model with dicranum
mixedmodel.full <- lmer(Shannon ~ Genus/Species + Compartment + 
                          SoilTemp + SoilpH + 
                          Genus:Compartment + 
                          Genus:SoilTemp +
                          Genus:SoilpH + 
                          Compartment:SoilTemp + 
                          Compartment:SoilpH + 
                          (1|Location/Microsite), data=div_plots); summary(mixedmodel.full)

sim_mod <- simulateResiduals(mixedmodel.full, refit=T, n=99)
plot(sim_mod)
Anova(mixedmodel.full)

# without dicranum
nodic_df <- div_plots %>% filter(Genus == "Polytrichum")

nodic_test <- nodic_df %>% mutate(across(c(is.character), factor))
mixedmodel.nodic <- lmer(Shannon ~ Species + Compartment + SoilTemp + SoilpH + 
                           Species:Compartment + 
                           Species:SoilTemp + 
                           Compartment:SoilTemp + Compartment:SoilpH +
                           (1|Microsite/Location), data=nodic_df); summary(mixedmodel.nodic)

sim_mod <- simulateResiduals(mixedmodel.nodic, refit=T, n=99)
plot(sim_mod)
Anova(mixedmodel.nodic)

plot_residuals(mixedmodel.full)

# r2 anova
lmfit <-  lm(model.response(model.frame(mixedmodel.nodic)) ~ fitted(mixedmodel.nodic))
summary(lmfit)$r.squared

# alpha-div graphs ####
(alphall <- 
div_plots %>% 
  ggplot(aes(x = Compartment, y = Shannon)) +
  geom_boxplot(aes(color = Compartment), lwd = 1, outlier.shape = NA) +
  geom_jitter(aes(fill = Compartment), position = position_jitter(0.3), shape = 21, size = 3, alpha = 0.6) +
  theme_bw() +
  ylab("Shannon Alpha-Diversity") + 
  scale_fill_manual('Section',values = compColours) +
  scale_colour_manual(values = compColoursContour) +
  xlab("SECTION") + ylab("SHANNON ALPHA-DIVERSITY") + 
  ggtitle("(a)") + theme_bw() +
  theme(title = element_text(size=16, color="black",face="bold"),
        axis.text.x  = element_text(size=14, color="black", angle=45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(size=14, color="black"),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size=14, color="black"),
        strip.text = element_text(size = 16, face = "italic"),
        panel.grid = element_blank(),
        legend.position = "none") + 
  guides(fill="none") + coord_fixed(3))


(alphacomp <- 
div_plots %>% 
    ggplot(aes(x = Compartment, y = Shannon)) +
    geom_boxplot(aes(color = Compartment), lwd = 1, outlier.shape = NA) +
    geom_jitter(aes(fill = Compartment), position = position_jitter(0.3), shape = 21, size = 4, alpha = 1) +
    facet_wrap("Host", ncol = 4, 
             labeller = label_wrap_gen(width=10)) +
  theme_bw() +
  ylab("Shannon Diversity") + 
  facet_wrap("Host", ncol = 4,
             labeller = labeller(label_wrap_gen(width=12),
                                 Host = mossNamesFormatted)) +
    scale_fill_manual('Section',values = compColours) +
    scale_colour_manual(values = compColoursContour) + theme_bw() +
  theme(title = element_text(size=16, color="black"),
        axis.text.x  = element_text(size=14, color="black", angle=45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(size=14, color="black"),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size=16, color="black"),
        strip.text = element_text(size = 18),
        panel.grid = element_blank(),
        legend.position = "none") + 
  guides(fill="none") ) #+ coord_fixed(3))
# use coord_fixed(3) for long

(alphaspec <- 
div_plots %>% 
  ggplot(aes(x = Host, y = Shannon)) +
  geom_boxplot(aes(color = Host), lwd = 1, outlier.shape = NA) +
  geom_jitter(aes(fill = Host), position = position_jitter(0.3), shape = 21, size = 3, alpha = 0.6) +
  facet_wrap("Compartment", ncol = 2, 
             labeller = label_wrap_gen(width=10)) +
  theme_bw() +
  ylab("Shannon Alpha-Diversity") + 
  scale_fill_manual(values=palette) +
  scale_x_discrete(labels = labelsItal) +
  scale_color_manual(values=palette) + 
  xlab("SECTION") + ylab("SHANNON ALPHA-DIVERSITY") + 
  ggtitle("(a)") + theme_bw() +
  theme(title = element_text(size=14, color="black",face="bold"),
        axis.text.x  = element_text(size=12, color="black", angle=45, vjust = 1, hjust = 1),
        axis.text.y  = element_text(size=12, color="black"),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size=12, color="black"),
        strip.text = element_text(size = 16),
        panel.grid = element_blank(),
        legend.position = "none") + 
  guides(fill="none") + coord_fixed(3)
)

########################
#### BETA DIVERSITY ####
########################

div_plots <- data.frame(sample_data(rare_MAGs),
                        Shannon=rare_rich$Shannon) %>% 
  mutate(Moss = Host) %>% 
  separate(Host, c("Genus", "Species"), "_") %>% 
  mutate(
    Genus = case_when(
      Genus == "P" ~ "Polytrichum",
      Genus == "D" ~ "Dicranum",
      TRUE ~ Genus
    )
  ) %>% dplyr::rename(Host = Moss)

nodic <- subset_samples(psMossMAGs, !Host=="D_undulatum")
nodic_meta <- as(sample_data(nodic), "data.frame")
deseq_MAGs <- phyloseq_to_deseq2(ps100, ~Compartment)

# Create function for geometric mean (taken from online website of phyloseq)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Estimate factor size
deseq_MAGs = estimateSizeFactors(deseq_MAGs, geoMeans = apply(counts(deseq_MAGs), 1, gm_mean))

# Blind version
vst <- DESeq2::varianceStabilizingTransformation(deseq_MAGs, blind=T)
vst.mat <- SummarizedExperiment::assay(vst) # Extract transformed asv table
t.vst.mat<-t(vst.mat)
t.vst.mat[which(t.vst.mat<0)]<-0
t.comm.vst <- vegdist(t.vst.mat, "bray")
t.pcoa.vst <- capscale(t.vst.mat~1, distance="bray")
t.pcoa.vst$CA$eig[1:3]/sum(t.pcoa.vst$CA$eig)
eig.test <- t.pcoa.vst$CA$eig

div_plots$PC1 <- scores(t.pcoa.vst)$sites[,1]
div_plots$PC2 <- scores(t.pcoa.vst)$sites[,2]

########################
# Full permanova ########
########################
perm <- how(nperm = 999)

setBlocks(perm) <- with(div_plots, Location)
permanova.full <-adonis2(t.vst.mat~ Genus/Species + Compartment + SoilTemp + SoilpH + 
                           Genus:Compartment + 
                           Genus:SoilTemp + 
                           Genus:SoilpH + 
                           Compartment:SoilTemp + 
                           Compartment:SoilpH + 
                           SoilTemp:SoilpH ,
                         data=div_plots,
                    permutations = perm, 	
                    by = "terms");permanova.full 

# no dicranum permanovas  ####
nodic <- subset_samples(ps100, !Host=="D_undulatum")
nodic_meta <- as(sample_data(nodic), "data.frame")
nodic_deseq <- phyloseq_to_deseq2(nodic, ~Compartment)

# Estimate factor size
nodic_deseq = estimateSizeFactors(nodic_deseq, geoMeans = apply(counts(nodic_deseq), 1, gm_mean))

# Blind version
nodic.vst <- DESeq2::varianceStabilizingTransformation(nodic_deseq, blind=T)
nodic.vst.mat <- SummarizedExperiment::assay(nodic.vst) # Extract transformed asv table
t.nodic.vst.mat<-t(nodic.vst.mat)
t.nodic.vst.mat[which(t.nodic.vst.mat<0)]<-0

t.nodic.comm.vst <- vegdist(t.nodic.vst.mat, "bray")
t.nodic.pcoa.vst <- capscale(t.nodic.vst.mat~1, distance="bray")
t.nodic.pcoa.vst$CA$eig[1:3]/sum(t.nodic.pcoa.vst$CA$eig)
nodic.eig.test <- t.nodic.pcoa.vst$CA$eig
nodic_df$MDS1<-scores(t.nodic.pcoa.vst)$sites[,1]
nodic_df$MDS2<-scores(t.nodic.pcoa.vst)$sites[,2]

setBlocks(perm) <- with(nodic_df, Location)
permanova.nodic <-adonis2(t.nodic.vst.mat~ Species + Compartment + 
                            SoilTemp + SoilpH + 
                            Species:Compartment + Species:SoilTemp + Species:SoilpH + 
                            Compartment:SoilTemp + Compartment:SoilpH + 
                            SoilTemp:SoilpH ,data=nodic_df,
                         permutations = perm, 	
                         by = "terms");permanova.nodic

# Beta div plots ####
# Ordination plot between HOST SPECIES 

(whole <- div_plots %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = Host, shape = Compartment), size=5, color = 'black') +
  stat_ellipse(level=0.9,  color='black', geom = "polygon", alpha = 0.18, 
               aes(fill = Host)) +
  scale_fill_manual(values=palette,
                    labels = mossSpecies) + 
  scale_shape_manual(values = c(22, 21)) + 
  scale_color_manual(values=palette) +
  guides(fill = guide_legend(override.aes = list(fill = palette,
                                                      color = palette)),
         shape = "none") +
  theme_bw() + 
  theme(plot.title = element_text(size = 18),
        title = element_text(size=14, color="black"),
        axis.text.x = element_text(size=14, color="black"), 
        axis.text.y  = element_text(size=14, color="black"),
        axis.title  = element_text(size=16, color="black"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text.align = 0,
        #legend.position = "none",
        legend.text = element_text(colour="black", size = 14)) )#+ coord_fixed())
# if keeping the shape legend, use this: shape = guide_legend(override.aes = list(shape = c(15,16)))

(facetg <- div_plots %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = Compartment, shape = Compartment), size=5, color = 'black') +
  stat_ellipse(level=0.9,  color='black', geom = "polygon", alpha = 0.18, 
               aes(fill = Compartment)) +
  scale_fill_manual('Section',values = compColours) +
  scale_colour_manual('Section', values = compColoursContour) +
  scale_shape_manual('Section', values = c(22,21)) +
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        axis.title.x = element_text(size=16, color="black"),
        axis.text.x = element_text(size=14, color="black"), 
        axis.text.y  = element_text(size=14, color="black"),
        axis.title.y  = element_blank(), #element_text(size=12, color="black"),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text.align = 0,
        legend.text = element_text(colour="black", size = 14),
        legend.position = "none",
        strip.text = element_text(size = 18)) + 
  facet_wrap("Host", labeller = labeller(label_wrap_gen(width=12),
                                         Host = mossNamesFormatted))) 

########################
