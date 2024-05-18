################################################################################
### AUTHOR : JONATHAN RONDEAU-LECLAIRE #########################################
### Differential abundance tests ###############################################
################################################################################

library(pacman)
p_load(phyloseq, ANCOMBC, betareg, tidyverse, magrittr, DESeq2, vegan, RColorBrewer, bestNormalize,
       vegan, ComplexHeatmap, colorRamp2, circlize, patchwork, MetBrewer)
source("myFunctions.R")
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
# write_rds(DA_pw_host_Order,"data/R_out/DA_pw_host_Order.RDS")


#########################################################################
#****************************** SANDBOX ****************************#####
#########################################################################

# global test 

df_fig_global <- 
  DA_global$res %>% 
  dplyr::select(taxon, contains("Host")) %>% 
  dplyr::left_join(DA_global$res_global %>% 
                     dplyr::transmute(taxon, 
                                      diff_Host = diff_abn,
                                      passed_ss = passed_ss)) %>% 
  dplyr::filter(diff_Host == 1) %>% 
  dplyr::mutate(color = ifelse(passed_ss == 1, "black", "red")) %>%
  dplyr::transmute(taxon,
                   color = color,
                   `P. commune` = round(lfc_HostP_commune, 2),
                   `P. piliferum` = round(lfc_HostP_piliferum, 2), 
                   `P. juniperum` = round(lfc_HostP_juniperinum, 2)) %>%
  tidyr::pivot_longer(cols = `P. commune`:`P. juniperum`, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon) %>% 
    dplyr::mutate(group = factor(group, 
                                 levels = c("P. commune",
                                            "P. piliferum",
                                            "P. juniperum")))

lo = floor(min(df_fig_global$value))
up = ceiling(max(df_fig_global$value))
mid = (lo + up)/2
df_fig_global %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes for globally significant taxa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(color = df_fig_global %>%
                                     dplyr::distinct(taxon, color) %>%
                                     .$color))

# Which MAGs are the most abundant across all samples?
# Which MAGs are cyanobacteria?

# top_Brown <- df_fig %>%
#   arrange(lfc_CompartmentGreen) %>%
#   slice_head(n = 4) %>%
#   transmute(taxon = taxon,
#             lfc = lfc_CompartmentGreen)
# 
# top_Green <- df_fig %>%
#   arrange(lfc_CompartmentGreen) %>%
#   slice_tail(n = 4) %>%
#   transmute(taxon = taxon,
#             lfc = lfc_CompartmentGreen)



waterfall_plot <- function(sp_name, df, title, caption) {
  df_fig <<- df$res %>% 
    dplyr::select(taxon, ends_with(sp_name)) %>% 
    dplyr::filter(!!sym(paste0("diff_",sp_name)) == 1 &
                    !!sym(paste0("passed_ss_",sp_name)) == TRUE) %>% 
    dplyr::arrange(desc(!!sym(paste0("lfc_",sp_name))) ) %>%
    dplyr::mutate(direct = factor(ifelse(!!sym(paste0("lfc_",sp_name)) > 0, "Positive LFC", "Negative LFC"),
                                  levels = c("Positive LFC", "Negative LFC")),
                  taxon = factor(taxon, levels = unique(taxon))) 
  
  df_fig %>%
    ggplot(aes(x = taxon, y = !!sym(paste0("lfc_",sp_name)), fill = direct)) + 
    geom_bar(stat = "identity", width = 1, color = "white",
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = !!sym(paste0("lfc_",sp_name)) - !!sym(paste0("se_",sp_name)), 
                      ymax = !!sym(paste0("lfc_",sp_name)) + !!sym(paste0("se_",sp_name))), 
                  width = 0.2, position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log fold change", 
         title = title, caption=caption) + 
    scale_fill_manual(values = c("Positive LFC" = "springgreen4", "Negative LFC" = "lightsalmon4"))+
    scale_color_discrete(name = NULL) +
    theme_minimal(base_size = 20) + 
    theme(panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 70, hjust = 1,
                                     color = df_fig$color),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")) +
    guides(fill = FALSE)
}

waterfall_plot("CompartmentGreen", DA_pairwise_comp, 
               "", #Differentially abundant species by compartment.
               "Significantly abundant at p<0.05 (ajdusted).
               Restricted to species with >10% relative abundance that passed the sensitivity analysis.")


