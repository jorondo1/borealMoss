library(pacman)
p_load(ANCOMBC, betareg, tidyverse, magrittr, DESeq2, vegan, RColorBrewer, bestNormalize,
       phyloseq, vegan, ComplexHeatmap, colorRamp2, circlize, patchwork)
source("myFunctions.R")
moss.ps <- readRDS("data/R_out/mossMAGs.RDS")

############################################################
### Ancom-BC differential abundance across compartment ######
############################################################

DA_pairwise_comp <- readRDS("data/R_out/DA_pairwise_comp.RDS")
DA_pairwise_comp <- ancombc2(
  data = moss.ps, 
  tax_level= "Species",
  p_adj_method="holm", 
  prv_cut = 0.10, 
  fix_formula="Host + Compartment", 
  group = "Host", # specify group if >=3 groups exist, allows structural zero detection 
  struc_zero = TRUE,
  pairwise = TRUE,
  alpha = 0.01,
  verbose = TRUE,
  n_cl = 10 # cores for parallel computing
)
# write_rds(DA_pairwise_comp,"data/R_out/DA_pairwise_comp.RDS")

# Process DAA output ; could potentially be merged/simplified with next command
sfx <- "CompartmentGreen"
speciesLFC <- DA_pairwise_comp$res %>% 
  dplyr::select(taxon, ends_with(sfx)) %>% 
  dplyr::filter(!!sym(paste0("diff_",sfx)) == 1 &
                !!sym(paste0("passed_ss_",sfx)) == TRUE) %>% 
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
            SD = se_CompartmentGreen) %>% 
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
DA_host_species <- readRDS("data/R_out/DA_host_species.RDS")
DA_host_species <- moss.ps %>% 
  ancombc2(tax_level= "Species", fix_formula="Host + Compartment", group = "Host", 
           struc_zero = TRUE, dunnet = TRUE, verbose = TRUE, n_cl = 10, 
           alpha = 0.01, )
# write_rds(DA_host_species,"data/R_out/DA_host_species.RDS")

# Keep only taxa for which at least one differential test is significant AND 
# passed the sensitivity analysis. Dynamically produce a list of conditions: 
parse_DAA_results <- function(DAA) {
  # Extract all column suffixes
  suffixes <- DAA$res_dunn %>% 
    dplyr::select(starts_with("diff_")) %>% colnames %>% 
    str_replace("diff_Host","")
  
  # Create a character vector of conditions
  conditions <- purrr::map_chr(suffixes, ~ paste0(
    "(`diff_Host", .x, "` == TRUE & `passed_ss_Host", .x, "` == TRUE)"
    )) %>% paste(collapse = " | ")
  
  # evaluate this condition in a filter argument:
  DAA$res_dunn %>%
    dplyr::select(-starts_with('W_'), -starts_with('p_')) %>% 
    filter(eval(parse(text = conditions))) %>% 
    # pivot to a long dataset, with one line per pairwise test per taxa
    pivot_longer(cols = -taxon, 
               names_to = c(".value", "Group"), 
               names_pattern = "(lfc|se|q|diff|passed_ss)_(.+)", 
               values_drop_na = TRUE) %>% 
    mutate(across(c(lfc,se), ~ case_when(diff==FALSE ~ 0,
                                         TRUE~.x)),
           textcolour = case_when(lfc==0 ~ "grey95", TRUE ~ "black"),
           Group = str_remove(Group, "Host")) %>% 
    # add taxonomy
    left_join(moss.ps@tax_table %>% as.data.frame %>% tibble,
              join_by("taxon" == "Species"))
}
hostDA <- parse_DAA_results(DA_host_species)
write_rds(hostDA, 'data/R_out/DA_host.RDS')


#########################################################################
# write_rds(DA_global, 'DA_global.RDS')

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


