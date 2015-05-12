library(seqinr)

seq_pos <- read.fasta("amyloid_pos.fasta")
seq_neg <- read.fasta("amyloid_neg.fasta")

max_len <- max(sapply(c(seq_pos, seq_neg), length))

all_seqs <- t(sapply(c(seq_pos, seq_neg), function(i)
  c(i, rep(NA, max_len - length(i)))))
