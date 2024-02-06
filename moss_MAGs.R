library(pacman)
p_load(tidyverse, rstatix)

# Import MAG statistics
read_tsv("data/genome.stats") %>% 
  # fix variables of interest
  transmute(`L50 (kb)` = ctg_L50/1000, 
            `N50 (kb)` = ctg_N50/1000,
            n_contigs = n_contigs,
            GC = gc_avg, 
            MBP = contig_bp/1000000, 
            MAG = str_extract(filename, "[^/]+$")) %>% 
  left_join(read_tsv("data/novel_quality_scores.txt",
                     col_names = c('MAG', 'comp', 'cont', 'QS')),
            by = 'MAG') %>% 
  mutate(MAG = str_remove(MAG, ".fa"),
         QS = QS/100) %>% 
  write_tsv("data/R_out/MAG_summary.tsv")

# Improvement in sample containment 
dbNames <- c("GTDB r207", "GTDB + MAGs", "Genbank + MAGs")
raw <- read_delim("data/cntm_sum.txt", col_names = c("Sample", "db", "cntm"))
cntm.df <- raw %>% 
  mutate(db = case_when(db == "custom" ~ dbNames[2],
                        db == "genbank" ~ dbNames[3],
                        db == "gtdb" ~ dbNames[1]),
         db = factor(db, levels = dbNames))

ggplot(cntm.df, aes(x = db, y = cntm, fill=db)) +
  geom_violin() + theme_minimal() + guides(fill = 'none') +
  labs(y = 'Sample k-mer containment', x = 'Database')

cntm.df %>% filter(db == dbNames[1]) %$% cntm %T>% hist %>% shapiro.test
cntm.df %>% filter(db == dbNames[2]) %$% cntm %T>% hist %>% shapiro.test
cntm.df %>% filter(db == dbNames[3]) %$% cntm %T>% hist %>% shapiro.test

raw_wide <- raw %>% pivot_wider(names_from = db, values_from = cntm) 
wilcox.test(raw_wide$custom, raw_wide$genbank) #NS