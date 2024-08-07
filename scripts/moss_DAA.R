################################################################################
### AUTHOR : JONATHAN RONDEAU-LECLAIRE #########################################
### Differential abundance tests ###############################################
################################################################################

library(pacman)
p_load(phyloseq, ANCOMBC, betareg, tidyverse, magrittr, DESeq2, vegan, RColorBrewer, bestNormalize,
       vegan, ComplexHeatmap, colorRamp2, circlize, patchwork, MetBrewer)
source("scripts/myFunctions.R")
moss.ps <- readRDS("data/R_out/mossMAGs.RDS")

############################################################
### Ancom-BC differential abundance across compartment ######
############################################################

DA_pairwise_comp <- readRDS("data/R_out/DA_pairwise_comp.RDS")
DA_pairwise_comp <- ancombc2(
  data = moss.ps, 
  tax_level= "Species",
  prv_cut = 0.10, 
  fix_formula="Host + Compartment + SoilpH + SoilTemp", 
  # rand_formula = '(1|Location)',
  group = "Host", # specify group if >=3 groups exist, allows structural zero detection 
  struc_zero = TRUE,
  pairwise = FALSE,
  alpha = 0.05,
  verbose = TRUE,
  n_cl = 10 # cores for parallel computing
)
# write_rds(DA_pairwise_comp,"data/R_out/DA_pairwise_comp.RDS")

# Process DAA output ; could potentially be merged/simplified with next command
sfx <- "CompartmentGreen"
speciesLFC <- DA_pairwise_comp$res %>% 
  dplyr::select(taxon, ends_with(sfx)) %>% 
  dplyr::filter(!!sym(paste0("diff_",sfx)) == 1 &
                !!sym(paste0("passed_ss_",sfx)) == TRUE &
                !!sym(paste0("q_",sfx)) < 0.01) %>% 
  dplyr::arrange(desc(!!sym(paste0("lfc_",sfx))) ) %>%
  dplyr::mutate(direct = factor(ifelse(!!sym(paste0("lfc_",sfx)) > 0, "Positive LFC", "Negative LFC"),
                                levels = c("Positive LFC", "Negative LFC")),
                taxon = factor(taxon, levels = unique(taxon)))

# Prep&export table for plots
speciesLFC_comp <- speciesLFC %>% 
  transmute(LFC = abs(!!sym(paste0("lfc_",sfx))), 
            Species = taxon, 
            compAss = case_when(!!sym(paste0("lfc_",sfx))>0 ~ "Green",
                                !!sym(paste0("lfc_",sfx))<0 ~ "Brown"),
            SE = se_CompartmentGreen) %>% 
  right_join(moss.ps %>% # identifier \ species association table
               tax_table %>% as.data.frame %>% 
               select(Species) %>% rownames_to_column("MAG"),
             by = 'Species') %>% 
  filter(!is.na(compAss))

# Export DA dataframe
write_rds(speciesLFC_comp, 'data/R_out/speciesLFC_comp.RDS')

################################################
######## PAIRWISE DUNN'S TEST ON HOST ###########
################################################

# Location correlates with environmental variables:
samdat <- moss.ps %>% sample_data %>% data.frame

extract_r2 <- function(aov) {
  temp_aov <- aov %>% summary
  ssq <- temp_aov[[1]]$`Sum Sq`
  pval <- temp_aov[[1]]$`Pr(>F)`[1]
  message(paste("R^2 =",round(ssq[1]/(ssq[1]+ssq[2]),2), ", p=", formatC(pval, format = "e", digits = 2)))
}
extract_r2(aov(SoilTemp ~ Location, data = samdat))
extract_r2(aov(SoilpH ~ Host, data = samdat))

chisq.test(samdat$Location, samdat$Host) 

### DA test : 
DA_host_species <- readRDS("data/R_out/DA_host_species.RDS")
DA_host_species <- moss.ps %>% 
  ancombc2(tax_level= "Species", 
           fix_formula="Host + Compartment", 
           group = "Host", 
           alpha = 0.05,
           struc_zero = TRUE, 
           dunnet = TRUE, 
           verbose = TRUE, 
           n_cl = 10)
# write_rds(DA_host_species,"data/R_out/DA_host_species.RDS")

hostDA <- parse_DAA_results(DA_host_species, 'dunn', 0.01, 
                            'Host', 'Species', moss.ps)
write_rds(hostDA, 'data/R_out/DA_host_results.RDS')

# order DAA across hosts
DA_pw_host_Order <- ancombc2(
  data = moss.ps, 
  tax_level= "Order",
  prv_cut = 0.10, 
  fix_formula="Host + Compartment", 
  group = "Host", 
  struc_zero = TRUE,
  pairwise = TRUE,
  alpha = 0.05,
  verbose = TRUE,
  n_cl = 10 # cores for parallel computing
) 
write_rds(DA_pw_host_Order,"data/R_out/DA_pw_host_Order.RDS")
