colNames<- c("staxids", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "score", "length", "pident", "nident", "mismatch", "qseqid", "qlen", "sseqid")
read_delim("data/S-22-POLJUN-B_paired_1aa.blastout", col_names=colNames) %>% View
