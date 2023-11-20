## Importing to Phyloseq  ####

t_comptype <- phyloseq(sample_data(comptype),
                       otu_table(t_asv, taxa_are_rows = TRUE),
                       tax_table(as.matrix(t_taxa))); t_comptype

t_comptype <- readRDS("ps_comptype.RDS") 

# 6801 taxa and 131 samples 
# Remember that this is after removing sequences with fewer than 40 sequences


# ALPHA-DIVERSITY #####

t_comptype.rare <- rarefy_even_depth(t_comptype, sample.size = 260000, replace = FALSE, verbose = TRUE, trimOTUs = TRUE, rngseed = 4466); t_comptype.rare

# Shannon diveristy 
t_comptype.div <- estimate_richness(t_comptype.rare, measure = c("Shannon")) 
shapiro.test(t_comptype.div$Shannon)
t_comptype.norm<-bestNormalize(t_comptype.div$Shannon) #orderNorm
shapiro.test(t_comptype.norm$x.t)
t_comptype.div$Shannon.norm<-t_comptype.norm$x.t

t_comptype.div <- cbind(sample_data(t_comptype), t_comptype.div)

# Plotting boxplots 
t_comptype.bp <-ggplot(t_comptype.div, aes(Compartment, Shannon.norm)) + theme_bw() +
  geom_boxplot(aes(fill = Compartment)) + facet_wrap("Species", ncol = 4,  labeller = label_wrap_gen(width=10)) + 
  ggtitle("TAXONOMIC DIVERSITY OF MOSS COMPARTMENTS") + ylab("SHANNON ALPHA-DIVERSITY") + 
  theme(axis.text.x  = element_text (size=12, color="black", angle=45, vjust = 1, hjust = 1),
        axis.text.y  = element_text (size=12, color="black"),
        axis.title.x  = element_text(size=12, color="black"),
        axis.title.y  = element_text(size=12, color="black"),
        strip.text = element_text(size = 14)) + 
  scale_fill_manual(values = c("tan4", "springgreen4")) + 
  guides(fill="none") + coord_fixed(2); t_comptype.bp

#
#
# BETA-DIVERSITY ####

t_comptype.p2d <- phyloseq_to_deseq2(t_comptype, ~Compartment)

# Create function for computation (taken from online website of phyloseq)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Estimate factor size
t_comptype.p2d = estimateSizeFactors(t_comptype.p2d, geoMeans = apply(counts(t_comptype.p2d), 1, gm_mean))

# Blind version
t.comptype.vst <- DESeq2::varianceStabilizingTransformation(t_comptype.p2d, blind=T)
t.comptype.vst.mat <- SummarizedExperiment::assay(t.comptype.vst) # Extract transformed asv table
t.comptype.t.vst.mat<-t(t.comptype.vst.mat)
t.comptype.t.vst.mat[which(t.comptype.t.vst.mat<0)]<-0

t.comptype.t.comm.vst <- vegdist(t.comptype.t.vst.mat, "jaccard")
t.comptype.t.pcoa.vst <-capscale(t.comptype.t.vst.mat~1,distance="jaccard")
t.comptype.t.comm.vst
t.comptype.t.pcoa.vst$CA$eig[1:3]/sum(t.comptype.t.pcoa.vst$CA$eig)
t.comptype.eig.test <- t.comptype.t.pcoa.vst$CA$eig
t.comptype.eig.test[1]/sum(abs(t.comptype.eig.test)) 
t.comptype.eig.test[2]/sum(abs(t.comptype.eig.test))
t.comptype.eig.test[3]/sum(abs(t.comptype.eig.test))

t_comptype.div$MDS1<-scores(t.comptype.t.pcoa.vst)$sites[,1]
t_comptype.div$MDS2<-scores(t.comptype.t.pcoa.vst)$sites[,2]

# Plotting 
t_comptype.pcoa <- ggplot(t_comptype.div, aes(MDS1, MDS2, colour = Compartment)) + 
  stat_ellipse(level=0.9, geom = "polygon", alpha = 0.18, aes(fill = Compartment)) +   
  geom_point(size = 5, aes(shape = Species)) + 
  ggtitle("TAXONOMIC BETA-DIVERSITY BETWEEN SITE TYPES") +
  theme_bw() +
  theme(plot.title = element_text(size = 18),
        legend.title = element_text(colour="black", size=16, face="bold"),
        legend.text = element_text(colour="black", size = 14)) + 
  guides(fill="none") + 
  scale_fill_manual(values = c("tan4", "springgreen3")) + 
  scale_color_manual(values = c("lightsalmon4", "springgreen4")) +  
  coord_fixed(); t_comptype.pcoa



