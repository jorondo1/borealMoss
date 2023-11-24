library(pacman)
p_load(ANCOMBC, betareg, tidyverse, magrittr, DESeq2, vegan, RColorBrewer, bestNormalize,
       phyloseq, vegan, ComplexHeatmap, colorRamp2, circlize)
source("myFunctions.R")
psMossMAGs <- readRDS("data/psMossMAGs.RDS")

# DA Ancom-BC differential abundance across compartment:
DA_pairwise_comp <- readRDS("data/DA_pairwise_comp")
DA_pairwise_comp <- ancombc2(
  data = psMossMAGs, 
  tax_level= "Species",
  p_adj_method="holm", 
  prv_cut = 0.10, 
  fix_formula="Host + Compartment", 
#  fix_formula="Host + Compartment + SoilpH + AmbientTemp + SoilTemp", 
  group = "Host", # specify group if >=3 groups exist, allows structural zero detection 
  struc_zero = TRUE,
  pairwise = TRUE,
  alpha = 0.05,
  verbose = TRUE,
  n_cl = 10 # 10 cores for parallel computing
); DA_pairwise_comp$res %>% colnames
write_rds(DA_pairwise_comp,"data/DA_pairwise_comp")

# Our waterfall plot function:
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
    geom_bar(stat = "identity", width = 0.7, color = "black", 
             position = position_dodge(width = 0.4)) +
    geom_errorbar(aes(ymin = !!sym(paste0("lfc_",sp_name)) - !!sym(paste0("se_",sp_name)), 
                      ymax = !!sym(paste0("lfc_",sp_name)) + !!sym(paste0("se_",sp_name))), 
                  width = 0.2, position = position_dodge(0.05), color = "black") + 
    labs(x = NULL, y = "Log fold change", 
         title = title, caption=caption) + 
    scale_fill_discrete(name = NULL) +
    scale_color_discrete(name = NULL) +
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 70, hjust = 1,
                                     color = df_fig$color),
          plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt")) +
    guides(fill = FALSE)
}

waterfall_plot("CompartmentGreen", DA_pairwise_comp, 
               "Differentially abundant species in green compartment",
               "Significantly abundant at p<0.01 (ajdusted).
               Restricted to species with >10% relative abundance that passed the sensitivity analysis.")

# Export DA dataframe
write_rds(df_fig, 'data/DA_results.RDS')

# Pathway groups we are interested in
pwGroups_interest <- c("Aromatics degradation", "Carbon fixation", 
                       "LPS metabolism", "Methane metabolism",
                       "Nitrogen metabolism", "Photosynthesis",
                       "Plant pathogenicity", "Sulfur metabolism",
                       "Symbiosis")

pwGroups_interest <- c("Carbon fixation", "Methane metabolism")

# Parse Module Completeness table
pwComp <- read_tsv("data/microbeannotator_out/metabolic_summary__module_completeness.tab") %>%
  dplyr::rename(pwGroup = `pathway group`) %>% 
  filter(pwGroup %in% pwGroups_interest) %>% 
  rename_with(~str_remove_all(.x, "\\.faa\\.ko")) %>% 
  rename_with(~simplify_name(.x)) %>% 
  mutate(pwName = paste0(gsub(" ","_", pwGroup), "_",name))
# Some have no .faa on NCBI, needs manual translation...

# taxID to Species conversion table
taxID_list <- psMossMAGs@tax_table %>% 
  data.frame %>% 
  select(Species) %>%
  rownames_to_column("taxID") %>% 
  mutate(taxID = simplify_name(taxID)) %>% 
  filter(taxID %in% colnames(pwComp))

# Prepare metadata for Heatmap
DAspec <- df_fig %>% 
  transmute(Species = taxon,
            Compartment = case_when(lfc_CompartmentGreen>0 ~ "Green",
                                    lfc_CompartmentGreen<0 ~ "Brown"),
            LFC = lfc_CompartmentGreen) %>% 
  filter(Species %in% taxID_list$Species) %>% 
  left_join(taxID_list, by = "Species")

# Prepare values matrix for Heatmap
mat <- pwComp %>% 
  dplyr::select(pwName, any_of(DAspec$taxID)) %>% 
  # Filter out rows where the max of all values (except pwName) is 30% or less
  filter(apply(select(., -pwName), 1, max) > 30) %>% 
  # Remove rows full of zeros
  filter(rowSums(select(., -pwName)) != 0) %>%
  column_to_rownames('pwName') %>% as.matrix 

# Create Heatmap!
ht <- Heatmap(mat,
        name = "Pathway Completeness across Compartments",
        col = colorRamp2(c(0, 50, 100), c("beige", "red", "darkorchid4")),
        column_split = DAspec$Compartment,
        heatmap_legend_param = list(
          title = "Completeness",
          title_position = "lefttop",
          legend_direction = "horizontal"
        ),
#        clustering_distance_rows = 'maximum',
        top_annotation = HeatmapAnnotation(df = data.frame(Compartment = DAspec$Compartment),
                                           col = list(Compartment = c("Green" = "springgreen3", 
                                                              "Brown" = "tan4")),
                                           which = "col",  # 'col' or 'row' based on the orientation
                                           show_legend = FALSE,
                                           show_annotation_name = FALSE
                                           ),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 12),
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        cluster_columns = TRUE)  # To keep the order of taxon as in DAspec
draw(ht, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#####################################################
### BETA REGRESSION between PW% and compartment #####
#####################################################

# create a transposed data set that will be used for hypothesis testing
pwComp_t <- pwComp %>% 
  select(-name,-pwGroup, -pwName) %>% 
  pivot_longer(cols = -module, 
               names_to = "taxID", 
               values_to = "value") %>% 
  pivot_wider(names_from = module, 
              values_from = value, id_cols = taxID) %>% 
  # add genome data we're interested in ...
  inner_join(DAspec, ., by = 'taxID') %>%
  # First, remove any row that's only zeros
  filter(rowSums(select_if(., is.numeric)) != 0) %>% 
  # transform because betareg doesn't allow 1 or 0 (see https://stats.stackexchange.com/questions/48028/beta-regression-of-proportion-data-including-1-and-0 )
  mutate(across(where(is.numeric), ~ ((.x / 100) * (length(.)-1) + 0.5 )/length(.))) 

# BETA REGRESSION
coefficients_list <- list()
for (i in colnames(pwComp_t[,5:dim(pwComp_t)[2]])) {
  print(i)
  model_call <- list(formula = as.formula(paste0(i,' ~ Compartment + LFC')), 
                     data = pwComp_t, 
                     link = 'logit')
    coefficients_list[[i]] <- do.call("betareg", model_call) %>% 
    summary %$% coefficients %>% 
    as.data.frame %>% 
    rownames_to_column("variable") %>% 
    filter(variable!="(Intercept)") %>% 
    mutate(pathway = i) %>% 
    dplyr::rename(p = mean.Pr...z..)
}
do.call(rbind, coefficients_list) %>% 
  filter(variable == "CompartmentGreen") %>% 
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ '***',
                         p_adj < 0.01 ~ '**',
                         p_adj < 0.05 ~ '*',
                         TRUE ~'')) %>% 
  select(variable, mean.Estimate, p_adj, p, sig, pathway) %>% 
  left_join(pwComp %>% select(module,name,pwGroup), 
            by = join_by(pathway == module)) %>% View


#########################################################################

# Global test for host
# DA_global <- readRDS("DA_global.RDS")

DA_global <- ancombc2(
  data = psMossMAGs, 
  tax_level= "Species",
  p_adj_method="holm", 
  prv_cut = 0.05, 
  fix_formula="Host + Compartment", 
  group = "Host", # specify group if >=3 groups exist, allows structural zero detection 
  struc_zero = TRUE,
  neg_lb = TRUE,
  pseudo_sens = TRUE,
  global = TRUE,
  pairwise = TRUE,
  alpha = 0.01,
  n_cl = 10,
#  dunnet = TRUE, 
  verbose = TRUE
)

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


