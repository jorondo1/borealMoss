library(pacman)
p_load(tidyverse, rstatix, scales)

# Import MAG statistics
read_tsv("data/genome.stats") %>% 
  # fix variables of interest
  transmute(`L50 (kb)` = ctg_L50/1000, 
            `N50 (kb)` = ctg_N50/1000,
            `Number of contigs`= n_contigs,
            GC = gc_avg, 
            MBP = contig_bp/1000000, 
            MAG = str_extract(filename, "[^/]+$")) %>% 
  left_join(read_tsv("data/novel_quality_scores.txt",
                     col_names = c('MAG', 'comp', 'cont', 'QS')),
            by = 'MAG') %>% 
  mutate(MAG = str_remove(MAG, ".fa"),
         QS = QS/100) %>% 
  write_tsv("data/R_out/MAG_summary.tsv")

MAG_summary <- moss.ps@tax_table %>% as.data.frame %>% 
  rownames_to_column("MAG") %>% as_tibble %>% 
  inner_join(read_tsv("data/R_out/MAG_summary.tsv"),
             by = 'MAG')

# MAG summary table
MAG_summary %>% 
  select(MAG, Class, Genus, `L50 (kb)`, n_contigs, GC, MBP, comp, cont, QS) %>% 
  arrange(Class, Genus) %>% 
  mutate(n_contigs = number(n_contigs, accuracy = 1)) %>% 
  mutate(across(where(is.numeric), ~number(.x, accuracy = 0.01))) %>% 
  write_csv("out/nMAG_statistics_full.csv")

# MAG summary statistics
MAG_summary %>% 
  summarise(across(where(is.numeric),
                   ~ paste0(round(mean(.),2), " ± ", round(sd(.),2)))) %>% View

# Eremio table
MAG_summary %>% 
  filter(Phylum == 'Eremiobacterota') %>% 
  select(MAG, Class, Genus, `L50 (kb)`, n_contigs, GC, MBP, comp, cont, QS) %>% 
  arrange(desc(QS)) %>% 
  mutate(comp = number(comp, accuracy = 0.1),
         n_contigs = number(n_contigs, accuracy = 1),
         across(where(is.numeric), ~number(.x, accuracy = 0.01))) %>% 
  write_csv("out/nMAG_stats_eremio.csv")

# Baltobacterales distribution 
balto <- moss.ps %>% psmelt %>% 
  group_by(Sample, Order, Compartment, Host) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), 
            .groups = 'drop') %>% 
  group_by(Sample) %>%
  mutate(relAb = Abundance/sum(Abundance)) %>% ungroup %>% 
  filter(Order == 'Baltobacterales') %>% group_by(Host)

balto %>% 
  summarise(mean_relab = mean(relAb),
            sd_relab = sd(relAb),
            n = n(),
            min = min(relAb))

balto %>% ggplot(aes(x = Host, y = relAb, fill = Host))  +
  geom_violin() + geom_jitter() + theme_light() +
  facet_wrap('Compartment', ncol = 2) +
  labs(y = "Relative abundance", x = '')

### Metabolism of Eremiobacterota for dall-e prompt:
full_join(
  read_delim("data/metabolic_summary__module_completeness_missing.tab"),
  read_delim("data/metabolic_summary__module_completeness.tab"),
  by = c('module', 'name', 'pathway group')) %>% 
  dplyr::rename(pwGroup = `pathway group`) %>% 
  rename_with(~str_remove_all(.x, "\\.faa\\.ko")) %>% 
  mutate(pwName = paste0(gsub(" ","_", pwGroup), "_",name),
         across(where(is.numeric), ~./100)) %>% 
  select(module, name, pwGroup, eremio$MAG) %>% 
  filter(rowSums(select(., where(is.numeric))) != 0) %$%
  pwGroup %>% unique %>% intersect(c("Aromatics degradation", "Carbon fixation", 
                                     "LPS metabolism", "Methane metabolism",
                                     "Nitrogen metabolism", "Photosynthesis",
                                     "Sulfur metabolism",
                                     "Symbiosis"))

