# ANCOMBC2 primer
# Compute Log-fold changes across 2+ groups of samples
library(pacman)
p_load(ANCOMBC, magrittr, dplyr, ggplot2, MetBrewer)

DA_Order <- ancombc2(
  data = phyloseq_object, 
  tax_level = "Order",
  prv_cut = 0.10, # Prevalence filtering
  fix_formula = "Treatment",
  # rand_effect = "Effet_aleatoire",
  group = "Treatment", # Groups to compare
  struc_zero = TRUE, # Structural zero detection 
  pairwise = TRUE, # Global test (all test types are FALSE by default)
  alpha = 0.05, 
  verbose = TRUE,
  n_cl = 10 # cores for parallel computing
)

# Parse the DAA results for a pairwise or dunnet test. Keep taxa with at 
# least one significant DA and having passed the sensitivity analysis. 
# !! Does NOT work with the global test, which has no LFC values.

parse_DAA <- function(DAA, # ANCOMBC output
                              test, # 'pair' or 'dunnet'
                              thr, # p-value threshold
                              gr, # group variable name (string)
                              taxRank) { # Taxonomic rank tested
  
  # define results list item name 
  listItem <- paste0('res_', test)
  
  # Extract suffixes for every relevant column
  suffixes <- DAA[[listItem]] %>% 
    dplyr::select(starts_with("diff_")) %>% colnames %>% 
    str_replace(paste0("diff_",gr),"")
  
  # Create a character vector of conditions
  conditions <- purrr::map_chr(suffixes, ~ paste0(
    "(`diff_", gr, .x, # passed sensitivity analysis:
    "` == TRUE & `passed_ss_", gr, .x, # adj. p value above certain threshold:
    "` == TRUE & `q_", gr, .x, "` < ", thr, ")" 
  )) %>% paste(collapse = " | ")
  
  DAA[[listItem]] %>%
    # drop irrelevant columns: 
    dplyr::select(-starts_with('W_'), -starts_with('p_')) %>% 
    # evaluate all these conditions through filter():
    filter(eval(parse(text = conditions))) %>% 
    # pivot to a long dataset, with one line per test per taxa
    pivot_longer(cols = -taxon, 
                 names_to = c(".value", "Group"), 
                 names_pattern = "(lfc|se|q|diff|passed_ss)_(.+)", 
                 values_drop_na = TRUE) %>% 
    mutate( # change LFC values to 0 if not significant (for the plot):
      across(c(lfc,se), ~ case_when(diff == FALSE ~ 0, TRUE ~ .x)),
      # plot colours:
      textcolour = case_when(lfc==0 ~ "white", TRUE ~ "black"),
      # format group names (needs improvement...)
      Group = str_remove(Group, gr)) %>% 
    # add taxonomy to the data (to allow colouring by higher rank)
    left_join(moss.ps@tax_table %>% as.data.frame %>% tibble %>% 
                dplyr::select(Domain:everything()[which(names(.) == taxRank)]) %>% 
                unique,
              by = c("taxon" = taxRank)) %>% 
    # remove all taxa for which there is no difference in any group:
    group_by(taxon) %>% 
    filter(!all(diff == FALSE)) %>% ungroup()
  # add sth similar to remove pairwise comparisons that have 0 DA taxa?
}

DA_order_results <- parse_DAA(DA_Order, 'pair', 0.01, 'Host', 'Order')

### Plot ! 
ggplot(data = DA_order_results,
       aes(x = Group, y = taxon, fill = lfc)) +
  geom_tile() +
  scale_fill_gradient2(low = met.brewer("Cassatt1")[1], 
                       mid = "white", 
                       high = met.brewer("Cassatt1")[8], 
                       midpoint = 0) +
  geom_text(aes(Group, taxon, label = round(lfc, 2), color = textcolour)) +
  scale_color_identity(guide = FALSE) + 
  theme_minimal() + 
  labs(x = '', y = '')




