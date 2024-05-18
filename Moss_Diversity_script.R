### AUTHOR : SARAH ISHAK ######################################################

library(pacman)
p_load(tidyverse, phyloseq, DESeq2, vegan, RColorBrewer, bestNormalize, wesanderson, sjPlot, car, lme4, emmeans, DHARMa, reshape2, SOfun, rstatix, ggpubr, circlize, picante, ggeffects,ComplexHeatmap)
#
#### Loading data ####
psMossMAGs <- readRDS("JO_psMossMAGs.RDS")
psMossPath <- as.data.frame(readRDS("psMossPathways_NEW.RDS"))

MAGs_melt <- psMossMAGs %>% psmelt %>% 
  select(OTU, Sample, Abundance, Compartment, Microsite, Host, Domain:Species)

asv <- as.data.frame(psMossMAGs@otu_table)

novelMAGs <- read.table("novel_MAGs.txt", sep=",")
MAGs_melt <- psMossMAGs@tax_table %>% data.frame %>% tibble::rownames_to_column(var = "V1")

# 
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
######################

# Shannon diversity calculated by rarefaction

rare_MAGs <- rarefy_even_depth(ps100,sample.size = 24900, replace = FALSE,verbose = TRUE,  trimOTUs = TRUE, rngseed = 4466)

rare_rich <- estimate_richness(rare_MAGs, measure = c("Shannon")) 

samsums_rare <- sample_sums(rare_MAGs)
samsums_rare[order(samsums_rare)]

shapiro.test(rare_rich$Shannon)
# the shapiro test came out normal for Shannon (0.024)...so we wouldn't need to normalise it...

div_plots <- data.frame(sample_data(rare_MAGs),
           Shannon=rare_rich$Shannon)

palette <- wes_palette("AsteroidCity1", n=4, type = "discrete")  

# modelling ####
div_plots <- div_plots %>% mutate(Moss = Host) %>% 
  separate(Host, c("Genus", "Species"), "_") %>% 
  mutate(
    Genus = case_when(
      Genus == "P" ~ "Polytrichum",
      Genus == "D" ~ "Dicranum",
      TRUE ~ Genus
    )
  ) %>% 
  dplyr::rename(Host = Moss)


div_plots <- div_plots %>% mutate(across(c(is.character), factor))


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

# no dicranum
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
r2.corr.mer <- function(m) {
  lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
  summary(lmfit)$r.squared
}

# Need to ask Isabelle 
r2.corr.mer(mixedmodel.nodic) 


# alpha-div graphs ####
(alphall <- 
div_plots %>% 
  ggplot(aes(x = Compartment, y = Shannon)) +
  geom_boxplot(aes(color = Compartment), lwd = 1, outlier.shape = NA) +
  geom_jitter(aes(fill = Compartment), position = position_jitter(0.3), shape = 21, size = 3, alpha = 0.6) +
  theme_bw() +
  ylab("Shannon Alpha-Diversity") + 
  scale_fill_manual('Section',values = c("darkgoldenrod4", "darkolivegreen3")) +
  scale_colour_manual(values = c("tan4", "springgreen4")) +
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
                                 Host = c("D_undulatum" = "D. undulatum",
                                             "P_commune" = "P. commune",
                                             "P_juniperinum" = "P. juniperinum",
                                             "P_piliferum" = "P. piliferum"))) +
    scale_fill_manual('Section',values = c("darkgoldenrod4", "darkolivegreen3")) +
    scale_colour_manual(values = c("tan4", "darkolivegreen4")) + theme_bw() +
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
  scale_x_discrete(labels = (c(expression(italic("D.  undulatum")),
                               expression(italic("P. commune")),
                               expression(italic("P.  juniperinum")),
                               expression(italic("P. piliferum"))))) +
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

div_plots <- data.frame(sample_data(rare_MAGs),
                        Shannon=rare_rich$Shannon)

div_plots <- div_plots %>% mutate(Moss = Host) %>% 
  separate(Host, c("Genus", "Species"), "_") %>% 
  mutate(
    Genus = case_when(
      Genus == "P" ~ "Polytrichum",
      Genus == "D" ~ "Dicranum",
      TRUE ~ Genus
    )
  ) %>% 
  dplyr::rename(Host = Moss)


nodic <- subset_samples(psMossMAGs, !Host=="D_undulatum")
nodic_meta <- as(sample_data(nodic), "data.frame")

#

deseq_MAGs <- phyloseq_to_deseq2(ps100, ~Compartment)

# Create function for computation (taken from online website of phyloseq)
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
t.comm.vst
t.pcoa.vst <- capscale(t.vst.mat~1, distance="bray")
t.pcoa.vst$CA$eig[1:3]/sum(t.pcoa.vst$CA$eig)
eig.test <- t.pcoa.vst$CA$eig
eig.test[1]/sum(abs(eig.test)) 
eig.test[2]/sum(abs(eig.test))
eig.test[3]/sum(abs(eig.test))

div_plots$PC1<-scores(t.pcoa.vst)$sites[,1]
div_plots$PC2<-scores(t.pcoa.vst)$sites[,2]

########################
# Full permanova ####

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

# Create function for computation (taken from online website of phyloseq)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Estimate factor size
nodic_deseq = estimateSizeFactors(nodic_deseq, geoMeans = apply(counts(nodic_deseq), 1, gm_mean))

# Blind version
nodic.vst <- DESeq2::varianceStabilizingTransformation(nodic_deseq, blind=T)
nodic.vst.mat <- SummarizedExperiment::assay(nodic.vst) # Extract transformed asv table
t.nodic.vst.mat<-t(nodic.vst.mat)
t.nodic.vst.mat[which(t.nodic.vst.mat<0)]<-0

t.nodic.comm.vst <- vegdist(t.nodic.vst.mat, "bray")
t.nodic.comm.vst
t.nodic.pcoa.vst <- capscale(t.nodic.vst.mat~1, distance="bray")
t.nodic.pcoa.vst$CA$eig[1:3]/sum(t.nodic.pcoa.vst$CA$eig)
nodic.eig.test <- t.nodic.pcoa.vst$CA$eig
nodic.eig.test[1]/sum(abs(nodic.eig.test)) 
nodic.eig.test[2]/sum(abs(nodic.eig.test))
nodic.eig.test[3]/sum(abs(nodic.eig.test))

nodic_df$MDS1<-scores(t.nodic.pcoa.vst)$sites[,1]
nodic_df$MDS2<-scores(t.nodic.pcoa.vst)$sites[,2]

perm <- how(nperm = 999)
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

palette <- wes_palette("AsteroidCity1", n=4, type = "discrete")  
 
(whole <- 
div_plots %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = Host, shape = Compartment), size=5, color = 'black') +
  stat_ellipse(level=0.9,  color='black', geom = "polygon", alpha = 0.18, 
               aes(fill = Host)) +
  scale_fill_manual(values=palette,
                    labels = (c("D. undulatum",
                                "P. commune",
                                "P. juniperinum",
                                "P. piliferum"))) + 
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



# if keeping the shape legend, use this 
# shape = guide_legend(override.aes = list(shape = c(15,16)))

(facetg <- 
div_plots %>% 
  ggplot(aes(PC1, PC2)) +
  geom_point(aes(fill = Compartment, shape = Compartment), size=5, color = 'black') +
  stat_ellipse(level=0.9,  color='black', geom = "polygon", alpha = 0.18, 
               aes(fill = Compartment)) +
  scale_fill_manual('Section',values = c("darkgoldenrod4", "darkolivegreen3")) +
  scale_colour_manual('Section', values = c("tan4", "springgreen3")) +
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
                                         Host = c("D_undulatum" = "D. undulatum",
                                                  "P_commune" = "P. commune",
                                                  "P_juniperinum" = "P. juniperinum",
                                                  "P_piliferum" = "P. piliferum")))) #+
  # coord_fixed())


########################

# Heatmap #### 
# Generate pw group index to split heatmap by row

pathways <- readRDS("psMossPathways_NEW.RDS")
DA_results <- readRDS("speciesLFC_comp.RDS")
psMossMAGs <- readRDS("JO_psMossMAGs.RDS")
DA_sub.tree <- readRDS("DA_sub.tree.RDS") #%>% dplyr::rename("taxon" = "Species")

#pathways <- pw 

# this is important, this is the order in which i want the pathway groups in
pathways2 <- pathways %>% mutate(across(`pathway group`, as_factor)) %>% 
  mutate(`pathway group` = fct_relevel(`pathway group`, c("Photosynthesis", 
                                                          "Nitrogen metabolism",
                                                          "Methane metabolism", 
                                                          "Prokaryotic carbon fixation"))) %>% 
  arrange(`pathway group`)
# check, and then 
pathways <- pathways2 


# transpose to get MAG name as rows, modules as columns :
pw_t <- pathways %>% dplyr::select(-module, -`pathway group`, -name) %>% t %>% 
  data.frame %>% setNames(pathways$module) %>% rownames_to_column(var = "id") 

species_list <- as.data.frame(psMossMAGs@tax_table) %>% select(Order, Species) %>% 
  rownames_to_column(var="id")
# you can do this to rename the new column directly select(Order, taxon = Species)

# we need to add the species list column because, this is what we will join by in DA_results

# let's select only modules that total more than zero 
# new <- pathways %>% mutate(sum = rowSums(select(., 3:1114))) %>% filter(sum > 0)

pathways_species <- left_join(pw_t, species_list, by = "id")
# looks good 


everything <- left_join(DA_results, pathways_species, by = "Species") %>% arrange(Order) %>% arrange(reverse(compAss)) 
# so freaking awesome bro 


DA_pathways <- everything %>% tibble::column_to_rownames(var = "Species") %>% 
  select(starts_with("M0")) 

pwgroup <- pathways %>% select(module, `pathway group`)

# load in pathway descriptions 
pw_desc <- pathways %>% select(module, name)


# let's select only columns that total more than zero 
new <- DA_pathways %>% select_if(colSums(.) > 0) %>% t %>% data.frame %>% 
  rownames_to_column(var = "module")

# make grouping variable thingy whatever 
new_pwgroup <- left_join(new, pwgroup) %>% select(module, `pathway group`)

pwgroup_vec <- new_pwgroup %>% column_to_rownames(var = "module") 
names2 <- rownames(pwgroup_vec) 
unlisted <- unlist(pwgroup_vec)
grouping = structure(unlisted, names = names2)

# here i'm trying to replace the module names with this, so it's more meaningful in the final output 
new_pwdesc <- left_join(new, pw_desc) %>% select(module, name)
new_pwdesc$name <- sub(" =>.*", "", new_pwdesc$name) # removing that => blah blah bit 

pwdesc_vec <- new_pwdesc %>% column_to_rownames(var = "module") 
pwdesc_names <- rownames(pwdesc_vec) 
unlisted_desc <- unlist(pwdesc_vec)
desc_labels = structure(unlisted_desc, names = pwdesc_names)

test <- new %>% tibble::column_to_rownames(var="module")

# let's see if i can't add an order category... 
order_vec <- everything %>% select(Species, Order) %>% column_to_rownames(var = "Species") %>% mutate(
  Order = case_when(
    Order == "Acidobacteriales" ~ "Terriglobales",
    TRUE ~ Order
  )
) %>% mutate()
nameso <- rownames(order_vec) 
unlistedo <- unlist(order_vec)
grouping_o = structure(unlisted, names = names2)


library(circlize)
col_fun = colorRamp2(c(100, 80, 60, 20, 0), c("#2c2d54", "#434475", "#6b6ca3", "#969bc7", "white"))



col_order = list(Order = c(
  "Acetobacterales" = "#D4A6C8FF", "Acidobacteriales" = "#A0CBE8FF", 
  "Actinomycetales" = "#F28E2BFF", "Armatimonadales" = "#FFBE7DFF",
  "Baltobacterales" = "#FF9D9AFF", "Bryobacterales" = "#FABFD2FF",
  "Burkholderiales" = "#9D7660FF", "Caulobacterales" = "#34b6c6",   
  "Chitinophagales" = "#E15749FF", "Chthoniobacterales" = "#D7B5A6FF", 
  "Cyanobacteriales" = "#59A14FFF", "Cytophagales" = "#8CD17DFF", 
  "Enterobacterales_A" = "#79706EFF", "Ktedonobacterales" = "#BAB0ACFF", 
  "Methylacidiphilales" = "#D37295FF", "Mycobacteriales" = "#961f1f", 
  "Rhizobiales" = "#499894FF", "Solirubrobacterales" = "#59386c",
  "Sphingobacteriales" = "#B07AA1FF", "Sphingomonadales" = "#4E79A7", 
  "Steroidobacterales" = "#132b69", "Streptosporangiales" = "#86BCB6FF", 
  "Xanthomonadales" = "#F1CE63FF", 
  "Pseudomonadales" = "hotpink", "Reyranellales" = "#B6992DFF"
  ))




Heatmap(as.matrix(test), col=col_fun,
        column_split=everything$compAss,
        row_split = grouping, 
        column_names_gp = gpar(fontsize = 24),
        column_title_gp = gpar(fontsize = 34, fontface = "bold"),
        column_names_side = "bottom",
        column_names_rot = 36,
        column_gap = unit(5, "mm"),
        show_column_dend = FALSE,
        cluster_columns = FALSE,
        row_names_max_width = max_text_width(desc_labels),
        show_row_dend = FALSE,
        row_labels = desc_labels,
        row_gap = unit(5, "mm"),
        row_names_gp = gpar(fontsize = 26),
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 32),
        cluster_rows = FALSE, 
        cluster_row_slices = FALSE, 
        bottom_annotation = HeatmapAnnotation(df = data.frame(Order = everything$Order),
                                              which = "col", 
                                              col = col_order,
                                              height = unit(8, "cm"),
                                              show_legend = FALSE, 
                                              show_annotation_name = FALSE,
                                              border=TRUE,
                                              simple_anno_size = unit(2, "cm")),
        top_annotation = HeatmapAnnotation(df = data.frame(Compartment = everything$compAss),
                                           col = list(Compartment = c("Green" = "darkolivegreen3",
                                                                      "Brown" = "darkgoldenrod4")),
                                           #which = "col",  # 'col' or 'row' based on the orientation
                                           show_legend = FALSE,
                                           show_annotation_name = FALSE,
                                           simple_anno_size = unit(2, "cm")),
        width = unit(90, "cm"),
        height = unit(55, "cm"),
        show_heatmap_legend = FALSE,
        border = TRUE
        ) # To keep the order of taxon as in DAspec)
# export at width = 5800, height = 2800


# 'darkgoldenrod4', 'darkolivegreen3'

## Legends ####
legend_list = list(Order = c(
  "Acetobacterales" = "#D4A6C8FF", "Terriglobales" = "#A0CBE8FF", "Burkholderiales" = "#9D7660FF",
  "Cyanobacteriales" = "#59A14FFF", "Reyranellales" = "#B6992DFF", "Rhizobiales" = "#499894FF",
  "Solirubrobacterales" = "#59386c", "Sphingomonadales" = "#4E79A7", "Steroidobacterales" = "#132b69",
  "Streptosporangiales" = "#86BCB6FF", "Xanthomonadales" = "#F1CE63FF"
))

colours <- unlist(legend_list)
taxa.names <- names(legend_list[[1]]) 

lgd2 = Legend(labels = taxa.names, title = "Order", 
              legend_gp = gpar(fill = colours))

lgd = Legend(col_fun = col_fun, title = "Pathway \nCompleteness (%)", border = "black")


draw(lgd, just = c("right"))
draw(lgd2, just = c("left"))


# make sure to clear plots first 
# wow this was sooooo fun thank you complexheatmaps for being easier than circle plots 


########################
