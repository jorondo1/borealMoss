library(pacman)
p_load(phyloseq, tidyverse, magrittr, ComplexHeatmap, colorRamp2, circlize,
       purrr, patchwork)
source("myFunctions.R")

#####################
##### METABOLISM #####
#####################

# Pathway groups we are interested in
pwGroups_interest <- c("Aromatics degradation", "Carbon fixation", 
                       "LPS metabolism", "Methane metabolism",
                       "Nitrogen metabolism", "Photosynthesis",
                       "Sulfur metabolism",
                       "Symbiosis")

pwGroups_interest <- c("Carbon fixation", "Nitrogen metabolism",
                       "Photosynthesis")

# Parse Module Completeness table
pwComp <- full_join(
  read_delim("data/metabolic_summary__module_completeness_missing.tab"),
  read_delim("data/metabolic_summary__module_completeness.tab"),
  by = c('module', 'name', 'pathway group')) %>% 
  dplyr::rename(pwGroup = `pathway group`) %>% 
  filter(pwGroup %in% pwGroups_interest) %>% 
  rename_with(~str_remove_all(.x, "\\.faa\\.ko")) %>% 
  rename_with(~simplify_name(.x)) %>% 
  mutate(pwName = paste0(gsub(" ","_", pwGroup), "_",name),
         across(where(is.numeric), ~./100))
# Some have no .faa on NCBI, needs manual translation...

# taxID to Species conversion table
taxID_list <- moss.ps@tax_table %>% 
  data.frame %>% 
  select(Species) %>%
  rownames_to_column("taxID") %>% 
  mutate(taxID = simplify_name(taxID)) %>% 
  filter(taxID %in% colnames(pwComp))

# Prepare metadata for Heatmap
DAspec <- readRDS('data/R_out/speciesLFC_comp.RDS') %>% 
  left_join(taxID_list, by = "Species")

# Prepare values matrix for Heatmap
pattern_regex <- gsub(" ", "_", pwGroups_interest) %>% 
  paste0("_") %>% 
  paste(collapse = "|")

prep_mat <- pwComp %>% 
  dplyr::select(pwName, any_of(DAspec$taxID)) %>% 
  # Filter out rows where the max of all values (except pwName) is 30% or less
  filter(apply(select(., -pwName), 1, max) > 0.30) %>% 
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
ht <- Heatmap(
  mat,
  name = "Pathway Completeness across Compartments",
  col = colorRamp2(c(0, 0.5, 1), c("beige", "red", "darkorchid4")),
  column_split = DAspec$compAss,
  row_split = groupIndex,
  heatmap_legend_param = list(
    title = "Completeness",
    title_position = "lefttop",
    legend_direction = "horizontal"),
  top_annotation = HeatmapAnnotation(
    df = data.frame(Compartment = DAspec$compAss),
    col = list(Compartment = c("Green" = compColours[2], 
                               "Brown" = compColours[1])),
    which = "col",  # 'col' or 'row' based on the orientation
    show_legend = FALSE,
    show_annotation_name = FALSE),
  row_names_gp = gpar(fontsize = 8),
  #row_names_rot=-45,
  column_names_gp = gpar(fontsize = 12),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  cluster_columns = TRUE)  # To keep the order of taxon as in DAspec

draw(ht, heatmap_legend_side = "bot", 
     annotation_legend_side = "bot",
     #  column_title = "Carbon fixation and Methane metabolism pathway completeness.",
     # column_title_gp=grid::gpar(fontsize=16)
)

#################################################
### Relationship between PW% and compartment #####
#################################################

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
  left_join(read_tsv("data/R_out/MAG_summary.tsv") %>% 
              select(MAG, comp), by = 'MAG') %>% 
  mutate(comp = case_when(is.na(comp) ~ 100, TRUE ~ comp)) %>% 
  # First, remove any row that's only zeros
  filter(rowSums(select_if(., is.numeric)) != 0) %>% 
  # LFC needs to be absolute value for stat test : 
  ##### UNLESS we use either compAss OR LFC ?!
#  mutate(LFC = abs(LFC)) ; ncol(pwComp_t) %>% 
  tibble

# remove columns where more than 90% of values are lower than 0.20
removeCols <- sapply( # apply to columns matching pattern name:
  pwComp_t[pwComp_t %>% names %>% grep("^M.{5}$", ., value = TRUE)], 
  # extract column names: 
  function(x) {sum(x <= 0.50) / length(x) >= 0.90}) %>% names(.)[.]

# Step 3: Remove identified columns from the dataframe
pwComp_t %<>%  select(-all_of(removeCols))
print(paste(ncol(pwComp_t), "columns left"))

# Regress !
coefficients_list <- list()
for (i in pwComp_t %>% names %>% grep("^M.{5}$", ., value = TRUE)) {
  #print(i)
  model_call <- list(formula = as.formula(paste0(i,' ~ compAss + comp')), 
                     data = pwComp_t, 
                     # use a quasi-binomial, rationale here 
                     family = quasibinomial(link = "logit")
  )
  coefficients_list[[i]] <- do.call("glm", model_call) %>% 
    summary %$% coefficients %>% 
    as.data.frame %>% 
    rownames_to_column("variable") %>% 
    filter(variable!="(Intercept)") %>% 
    mutate(pathway = i) %>% 
    dplyr::rename(p = `Pr(>|t|)`)
}
results <- do.call(rbind, coefficients_list) %>% 
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ '***',
                         p_adj < 0.01 ~ '**',
                         p_adj < 0.05 ~ '*',
                         TRUE ~'')) %>% 
  select(variable, Estimate, p_adj, p, sig, pathway) %>% 
  left_join(pwComp %>% select(module,name,pwGroup), 
            by = join_by(pathway == module))

results %>% filter(p_adj<0.05) %>% select(-p, -sig) %>%  View

# transpose to get MAG name as rows, modules as columns :
# pw_t <- pw %>% dplyr::select(-module) %>% t %>% 
#   data.frame %>% setNames(pw$module)


glm(formula = as.formula('M00374 ~ comp'), 
    data = pwComp_t,
    family = quasibinomial(link = "logit")) %>% summary



# DA results only list the Species, not the shortened MAG name
# Pathways matrix only has shortened MAG name as rows + modules as columns
# I extracted the Species for these MAGs from the taxonomy list in JO_psMossMAGs.RDS file
# But the MAG names written in the taxonomy file are written in the long form, so I had to shorten the names with one of the functions you wrote
# I think I just x_joined (can't remember if it was left or inner) to match the Species to the pathway matrix names then left_joined that to DA results to get a matrix of all the DA output + modules
