library(pacman)
p_load(ANCOMBC, betareg, tidyverse, magrittr, DESeq2, vegan, RColorBrewer, bestNormalize,
       phyloseq, vegan, ComplexHeatmap, colorRamp2, circlize, patchwork)
source("myFunctions.R")
moss.ps <- readRDS("data/psMossMAGs.RDS")

############################################################
### Ancom-BC differential abundance across compartment ######
############################################################

DA_pairwise_comp <- readRDS("data/DA_pairwise_comp")
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

# Waterfall plot
sfx <- "CompartmentGreen"
DA_comp.df <- DA_pairwise_comp$res %>% 
      dplyr::select(taxon, ends_with(sfx)) %>% 
      dplyr::filter(!!sym(paste0("diff_",sfx)) == 1 &
                      !!sym(paste0("passed_ss_",sfx)) == TRUE) %>% 
      dplyr::arrange(desc(!!sym(paste0("lfc_",sfx))) ) %>%
      dplyr::mutate(direct = factor(ifelse(!!sym(paste0("lfc_",sfx)) > 0, "Positive LFC", "Negative LFC"),
                                    levels = c("Positive LFC", "Negative LFC")),
                    taxon = factor(taxon, levels = unique(taxon)))

# Export DA dataframe
write_rds(DA_comp.df, 'data/DA_results.RDS')

################################################
######## PAIRWISE DUNN'S TEST ON HOST ###########
################################################
DA_host_species <- readRDS("data/DA_host_species.RDS")
DA_host_species <- moss.ps %>% 
  ancombc2(tax_level= "Species", fix_formula="Host + Compartment", group = "Host", 
           struc_zero = TRUE, dunnet = TRUE, verbose = TRUE, n_cl = 10, 
           alpha = 0.01, )
# write_rds(DA_host_species,"data/DA_host_species.RDS")

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

#####################
##### METABOLISM #####
#####################
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
taxID_list <- moss.ps@tax_table %>% 
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
pattern_regex <- gsub(" ", "_", pwGroups_interest) %>% 
  paste0("_") %>% 
  paste(collapse = "|")

prep_mat <- pwComp %>% 
  dplyr::select(pwName, any_of(DAspec$taxID)) %>% 
  # Filter out rows where the max of all values (except pwName) is 30% or less
  filter(apply(select(., -pwName), 1, max) > 30) %>% 
  # Remove rows full of zeros
  filter(rowSums(select(., -pwName)) != 0)

# Generate pw group index to split heatmap by row
groupIndex <- str_extract(prep_mat$pwName, pattern_regex) %>% 
  str_replace("_", " ") %>% str_replace("_","") %>% as.vector

# remove pw group from pathway name and generate final matrix
mat <- prep_mat %>%
  mutate(pwName = str_remove(pwName, pattern_regex)) %>% 
  column_to_rownames('pwName') %>% as.matrix

# Create Heatmap!
ht <- Heatmap(mat,
        name = "Pathway Completeness across Compartments",
        col = colorRamp2(c(0, 50, 100), c("beige", "red", "darkorchid4")),
        column_split = DAspec$Compartment,
        row_split = groupIndex,
        heatmap_legend_param = list(
          title = "Completeness",
          title_position = "lefttop",
          legend_direction = "horizontal"
        ),
        top_annotation = HeatmapAnnotation(df = data.frame(Compartment = DAspec$Compartment),
                                           col = list(Compartment = c("Green" = compColours[2], 
                                                              "Brown" = compColours[1])),
                                           which = "col",  # 'col' or 'row' based on the orientation
                                           show_legend = FALSE,
                                           show_annotation_name = FALSE
                                           ),
        row_names_gp = gpar(fontsize = 8),
        row_names_rot=-45,
        column_names_gp = gpar(fontsize = 12),
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        cluster_columns = TRUE)  # To keep the order of taxon as in DAspec

draw(ht, heatmap_legend_side = "bot", 
     annotation_legend_side = "bot",
   #  column_title = "Carbon fixation and Methane metabolism pathway completeness.",
    # column_title_gp=grid::gpar(fontsize=16)
   )

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
                     data = pwComp_t#, link = 'logit'
                     )
    coefficients_list[[i]] <- do.call("lm", model_call) %>% 
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


