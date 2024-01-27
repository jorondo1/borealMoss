### Moss compartment colours
compColours <- c('darkgoldenrod4', 'darkolivegreen3')

### Remove .__ from taxa names
rename.fun <- function(x) {str_remove(x,".__")}

### Parse Sourmash output
parse_SM <- function(gather_files) {
  Sys.glob(gather_files) %>% 
    map_dfr(read_csv, col_types="ddddddddcccddddcccd") %>%
    mutate(
      uniqueK = (unique_intersect_bp/scaled)*average_abund,
      genome = str_replace_all(name, c(" .*" = "", ".fa.sig" = "")), # remove _1 from MAG identifiers
      run = str_replace(query_name, "_clean", ""), 
      .keep="unused") %>% 
    select(run, uniqueK, genome) %>% 
    pivot_wider(names_from = run,
                values_from = uniqueK) %>% 
    replace(is.na(.), 0) %>% 
    arrange(genome) %>% 
    column_to_rownames("genome") %>% 
    round(digits=0)
}


# Renaming MAGs in a more compact way
simplify_name <- function(name) {
  pattern = "^([A-Za-z]{2})[^_]*_([^.]*)\\.bin\\.(\\d+)$"
  replacement = "\\1\\2\\3.MAG"
  str_replace_all(name, pattern, replacement)
}

# Make a phyloseq object 
makePhyloSeq <- function(abund, meta, tax, tree = FALSE) {
  ps <- phyloseq(
    otu_table(abund, taxa_are_rows = T),
    sample_data(meta),
    tax_table((tax %>% 
                 filter(genome %in% rownames(abund)) %>% 
                 column_to_rownames("genome")) %>% as.matrix)
    )
  if (tree!=FALSE) {ps <- read.tree(tree) %>% phyloseq %>% merge_phyloseq(ps)}
  return(ps)
}

# Prune
prunePS <- function(ps, keep) {
  ps %>% taxa_names %>% .[(. %in% keep)] %>% prune_taxa(ps)
}

trans.fun <- function(ps, method, keep) {
  ps %>% tax_filter(min_prevalence=0.05) %>%
    purrr::when(method == 'none' ~ ., # if satisfied, go to next pipe, otherwise
                ~ microbiome::transform(.,transform=method)) %>% # that function we defined above
    prunePS(.,keep)
}
