library(pacman)
p_load(tidyverse, magrittr, taxize, doParallel, foreach)
numCores <- detectCores()
registerDoParallel(numCores/2)
source("myFunctions.R")

colNames<- c("staxids", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", "qseqid", "qlen", "sseqid")
blastout <- read_delim("data/S-3-POLCOM-B_paired_1hd.blastout", col_names=colNames)

taxids <- blastout %$% staxids %>% unique

# don't use parallel processing, it generates error because it queries ncbi too much
# results <- foreach(id = taxids, .packages = "taxize") %dopar% {
#   tryCatch({
#     classification(id, db = "ncbi")
#   }, error = function(e) NA)
#                    }

results2 <- lapply(taxids, function(id) {
  tryCatch({
    classification(id, db = "ncbi")
  }, error = function(e) NA)
})

# Save that :
write_rds(results, "data/taxid_results")

taxLvls <- c("superkingdom", "phylum","class","order","family","genus","species")

taxonomy <- data.frame(matrix(ncol = length(taxLvls), nrow = 0))
colnames(taxonomy) <- taxLvls
for (i in 1:length(results)) {
  #df$taxid <- names(results[[i]][1])  # Add taxid as a new column
  taxonomy <- results[[i]][1][[1]] %>% 
    dplyr::filter(rank %in% taxLvls) %$% name %>% 
    c(names(results[[i]][1]),.) %>% 
  rbind(taxonomy)
}

