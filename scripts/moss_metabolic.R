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
  read_delim("data/metabolic_summary__module_completeness_missing.tab"),
  read_delim("data/metabolic_summary__module_completeness.tab"),
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
  mutate(across(where(is.numeric),  ~ round(.x, 3))) %>% 
  write_csv("out/metabolic_test_full.csv")
