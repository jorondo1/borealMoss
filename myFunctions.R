### Italicize species names
labelsItal <- c(expression(italic("D. undulatum")),
                expression(italic("P. commune")),
                expression(italic("P. juniperinum")),
                expression(italic("P. piliferum")))

### Order colours
col_order = c(
  "Acetobacterales" = "#4E79A7", "Acidobacteriales" = "#A0CBE8FF",
  "Actinomycetales" = "#F28E2BFF", "Armatimonadales" = "#FFBE7DFF",
  "Baltobacterales" = "#86BCB6FF", "Bryobacterales" = "#F1CE63FF",
  "Burkholderiales" = "#9D7660FF", "Caulobacterales" = "#34B6C6",
  "Chitinophagales" = "#E15749FF", "Chthoniobacterales" = "#D7B5A6FF",
  "Cyanobacteriales" = "#59A14FFF", "Cytophagales" = "#8CD17DFF",
  "Enterobacterales_A" = "#79706EFF", "Ktedonobacterales" = "#BAB0ACFF",
  "Methylacidiphilales" = "#D37295FF", "Mycobacteriales" = "#FABFD2FF",
  "Rhizobiales" = "#499894FF", "Solirubrobacterales" = "#59386C",
  "Sphingobacteriales" = "#B07AA1FF", "Sphingomonadales" = "#D4A6C8FF",
  "Steroidobacterales" = "#132B69", "Streptosporangiales" = "#FF9D9AFF",
  "Xanthomonadales" = "#961F1F",
  "Pseudomonadales" = "hotpink", "Reyranellales" = "#B6992DFF"
)

#####################
### Preprocessing ####
#####################

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

# Some taxonomic levels have redundancies because higher levels use alternative names,
# herego there can be a duplicate Order whose Class or Phylum is different. Some 
# can be heterotypic synonyms, others outdated taxonomic names.

listDupl <- function(tax, level) {
  level_sym <- rlang::sym(level)
  
  subset <- tax %>% dplyr::select(Domain:!!level_sym) %>% 
    dplyr::group_by(!!level_sym) %>%
    unique
  
  dupList <- subset %>% 
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n > 1) %>%
    dplyr::pull(!!level_sym)
  
  # Relative abundance from melted ps object 
  relab.fun <- function(df) {
    summarise(df, Abundance = sum(Abundance, na.rm = TRUE), 
              .groups = 'drop') %>% ungroup %>% 
      group_by(Sample) %>% 
      mutate(relAb = Abundance/sum(Abundance)) %>% ungroup
  }
  
  subset %>% arrange(!!level_sym) %>% 
    filter(!!level_sym %in% dupList) %>% print(n = 1000)
}

##################
### DIVERRSITY ####
##################

# Variance-stabilizing transformation
vst.fun <- function(ps, var) {
  phyloseq_to_deseq2(
    ps, as.formula(paste0('~', var))) %>% # DESeq2 object
    estimateSizeFactors(., geoMeans = apply(
      counts(.), 1, function(x) exp(sum(log(x[x>0]))/length(x)))) %>% 
    DESeq2::varianceStabilizingTransformation(blind=T) %>% # VST
    SummarizedExperiment::assay(.) %>% t %>% 
    { .[. < 0] <- 0; . } # replace negatives by zeros
}

# PCOA
pcoa.fun <- function(ps, var, vst.mx, dist) {
  # compute distance matrix
  dist.mx <- vegan::vegdist(vst.mx, distance = dist)
  # compute PCoA
  PCoA <- capscale(dist.mx~1, distance = dist)
  # Print first 3 coordinates
  eig <- round(PCoA$CA$eig[1:3]/sum(PCoA$CA$eig),2)
  message("First 3 principal coordinates :")
  message(paste(eig[1], ',', eig[2], ',', eig[3]))
  # create output list
  out <- data.frame(sample_data(ps))
  out$PCoA1 <- scores(PCoA)$sites[,1]
  out$PCoA2 <- scores(PCoA)$sites[,2]
  list(metadata = out, eig = PCoA$CA$eig, dissMx = vst.mx)
}

################
### PLOTTING ####
################

# Compute top taxa, aggregate residual in others
topTaxa <- function(psmelt, comp, taxLvl, topN) {
  psmelt %>% 
    filter(Compartment == comp) %>% # subset compartment
    group_by(!!sym(taxLvl)) %>% # group by tax level
    summarise(relAb = mean(relAb)) %>% # find top abundant taxa
    arrange(desc(relAb)) %>% 
    mutate(aggTaxo = as.factor(case_when( # aggTaxo will become the plot legend
      row_number() <= topN ~ !!sym(taxLvl), #+++ We'll need to manually order the species!
      row_number() > topN ~ 'Others'))) # +1 to include the Others section!
}

# Abundance summary dataset by compartment using top taxa
df_comm <- function(psmelt, comp, taxLvl, topTaxa) {
  # use topTaxa to join aggTaxo variable
  left_join(psmelt, topTaxa, by=taxLvl) %>% 
    dplyr::filter(Compartment == comp) %>% # subset
    aggregate(Abundance~aggTaxo+Host+Compartment, # Abundance is aggregated...
              data=., FUN = sum) %>% # ...sums "Others" taxa by specified variables
    mutate(aggTaxo = factor(aggTaxo,# reorder 
                            levels = topTaxaLvls)) 
}

# Extract tree labels in the order they appear :
pullTreeLabels <- function(plot, rank) {
  plot$data %>% filter(!is.na(label)) %>% 
    arrange(y) %>% select(!!sym(rank)) %>% 
    # dplyr way of collapsing a tibble into a vector is pull() :
    pull %>% unique
}
