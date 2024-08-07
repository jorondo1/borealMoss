################################################################################
### AUTHOR : JONATHAN RONDEAU-LECLAIRE #########################################
### Metabolism tests ###########################################################
################################################################################

library(pacman)
p_load(phyloseq, tidyverse, magrittr, ComplexHeatmap, colorRamp2, circlize,
       purrr, patchwork, kableExtra)
source("scripts/myFunctions.R")

##################################################################
##### PATHWAY COMPLETENESS AND COMPARTMENT-ASSOCIATED SPECIES #####
##################################################################

# Pathway groups we are interested in
pwGroups_interest <- c("Nitrogen metabolism", 
                       "Photosynthesis", 
                       "Methane metabolism", 
                       "Carbon fixation")

pwModules_interest <- c("M00173", "M00376", "M00375", "M00374",
                        "M00377", "M00579", "M00260")

# Parse Module Completeness table
pwComp <- full_join(
  read_delim("data/MAG_analysis/Annotation/microbeannotator_out/metabolic_summary__module_completeness_missing.tab"),
  read_delim("data/MAG_analysis/Annotation/microbeannotator_out/metabolic_summary__module_completeness.tab"),
  by = c('module', 'name', 'pathway group')) %>% 
  dplyr::rename(pwGroup = `pathway group`) %>% 
  # Keep modules of interest: 
  filter(pwGroup %in% pwGroups_interest | module %in% pwModules_interest) %>% 
  rename_with(~str_remove_all(.x, "\\.faa\\.ko")) %>% 
  rename_with(~simplify_name(.x)) %>% 
  mutate(pwName = paste0(gsub(" ","_", pwGroup), "_",name),
         across(where(is.numeric), ~./100))
# Some have no .faa on NCBI, needs manual translation...

#################################################
### Relationship between PW% and compartment #####
#################################################
# taxID to Species conversion table
taxID_list <- moss.ps@tax_table %>% 
  data.frame %>% 
  select(Species) %>%
  rownames_to_column("taxID") %>% 
  mutate(taxID = simplify_name(taxID)) %>% 
  filter(taxID %in% colnames(pwComp))

DAspec <- readRDS('data/R_out/speciesLFC_comp.RDS') %>% 
  left_join(taxID_list, by = "Species")

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
  left_join(read_tsv("data/MAG_analysis/MAG_summary.tsv") %>% 
              select(MAG, comp), by = 'MAG') %>% 
  mutate(comp = case_when(is.na(comp) ~ 100, TRUE ~ comp)) %>% 
  # First, remove any row that's only zeros
  filter(rowSums(select_if(., is.numeric)) != 0) %>% 
  # LFC needs to be absolute value for stat test : 
  ##### UNLESS we use either compAss OR LFC ?!
#  mutate(LFC = abs(LFC)) ; ncol(pwComp_t) %>% 
  tibble

# remove columns where more than 80% of values are lower than 0.20
removeCols <- sapply( # apply to columns matching pattern name:
  pwComp_t[pwComp_t %>% names %>% grep("^M.{5}$", ., value = TRUE)], 
  # extract column names: 
  function(x) {sum(x <= 0.20) / length(x) >= 0.80}) %>% names(.)[.]

# Remove identified columns from the dataframe
pwComp_t %<>% select(-all_of(removeCols))
print(paste(ncol(pwComp_t)-6, "columns left"))

# Overdispersion?
pwComp_t %>% dplyr::select(starts_with('M0')) %>% 
  pivot_longer(cols = starts_with('M'),
               values_to = 'Completeness',
               names_to = 'Module') %>% 
  group_by(Module) %>% 
  summarise(mean = mean(Completeness),
            sd = sd(Completeness)) %>% 
  mutate(diff = mean-sd) %>% 
  ggplot(aes(x = Module)) +
  geom_point(aes(y = mean), colour = 'blue') +
  geom_point(aes(y = sd), colour = 'red')

# Regress !
coefficients_list <- list()
for (i in pwComp_t %>% names %>% grep("^M.{5}$", ., value = TRUE)) {
  #print(i)
  model_call <- list(formula = as.formula(paste0(i,' ~ compAss + comp')), 
                     data = pwComp_t, 
                     # use a quasi-binomial
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
test_results <- do.call(rbind, coefficients_list) %>% 
  mutate(p_adj = p.adjust(p, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ '***',
                         p_adj < 0.01 ~ '**',
                         p_adj < 0.05 ~ '*',
                         TRUE ~'')) %>% 
  select(variable, Estimate, p_adj, p, sig, pathway) %>% 
  left_join(pwComp %>% select(module,name,pwGroup), 
            by = join_by(pathway == module))

test_results %>% filter(p_adj < 0.05 & variable == 'compAssGreen') %>% 
  select(-sig) %>% 
  write_csv("out/metabolic_test_sig.csv")

test_results %>% 
  filter(variable == 'compAssGreen') %>%  
  dplyr::select(pathway, name, pwGroup, Estimate, p_adj, sig) %>% 
  mutate(across(where(is.numeric),  ~ round(.x, 3))) %>% View
  write_csv("out/metabolic_test_full.csv")


#############
# Heatmap #### 
#############=
# Author : Sarah Ishak
# Second half of this code should be moved to moss_plots.R

mossPW <- as.data.frame(readRDS("data/R_out/pathway_completeness_table.RDS")) 
moss.ps <- readRDS("data/R_out/mossMAGs.RDS")
DA_results <- readRDS("data/R_out/speciesLFC_comp.RDS")

# this is important, this is the order in which i want the pathway groups in
pathways <- mossPW %>% 
  mutate(across(`pathway group`, as_factor)) %>% 
  mutate(`pathway group` = fct_relevel(`pathway group`, c("Photosynthesis", 
                                                          "Nitrogen metabolism",
                                                          "Methane metabolism", 
                                                          "Prokaryotic carbon fixation"))) %>% 
  arrange(`pathway group`)

# transpose to get MAG name as rows, modules as columns :
pw_t <- pathways %>% dplyr::select(-module, -`pathway group`, -name) %>% t %>% 
  data.frame %>% setNames(pathways$module) %>% rownames_to_column(var = "id") 

species_list <- as.data.frame(moss.ps@tax_table) %>% select(Order, Species) %>% 
  rownames_to_column(var="id")
# you can do this to rename the new column directly select(Order, taxon = Species)

# we need to add the species list column because, this is what we will join by in DA_results

# let's select only modules that total more than zero 
# new <- pathways %>% mutate(sum = rowSums(select(., 3:1114))) %>% filter(sum > 0)

pathways_species <- left_join(pw_t, species_list, by = "id")
# looks good 

everything <- left_join(DA_results, pathways_species, by = "Species") %>% 
  arrange(Order) %>% arrange(reverse(compAss)) 
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
order_vec <- everything %>% select(Species, Order) %>% 
  column_to_rownames(var = "Species") %>% 
  mutate(
    Order = case_when(
      Order == "Acidobacteriales" ~ "Terriglobales",
      TRUE ~ Order))

nameso <- rownames(order_vec) 
unlistedo <- unlist(order_vec)
grouping_o = structure(unlisted, names = names2)

library(circlize)
col_fun = colorRamp2(c(100, 80, 60, 20, 0), c("#2c2d54", "#434475", "#6b6ca3", "#969bc7", "white"))

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

