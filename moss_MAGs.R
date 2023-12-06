library(pacman)
p_load(phyloseq, tidyverse, magrittr, kable, kableExtra)
source("myFunctions.R")

MAG_summary <- inner_join(
  # Final set of novel taxa found in at least one sample:
  psMossMAGs@tax_table %>% data.frame %>% 
    rownames_to_column("MAG"),
  # MAG metrics
  read_delim("data/novel_quality_scores.txt",
             col_names = c('MAG', 'Completeness', 'Contamination', 'QualityScore')) %>% 
            mutate(MAG = str_remove(MAG, '.fa')),
          by = "MAG") %>% selecet(-Species)

MAG_summary %>% write_tsv("out/MAG_summary.txt")

