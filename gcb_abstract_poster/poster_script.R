library(seqinr)

seq_pos <- read.fasta("gcb_abstract_poster/amyloid_pos_full.fasta", seqtype = "AA")
seq_neg <- read.fasta("gcb_abstract_poster/amyloid_neg_full.fasta", seqtype = "AA")

max_len <- max(sapply(c(seq_pos, seq_neg), length))

all_seqs <- t(sapply(c(seq_pos, seq_neg), function(i)
  c(i, rep(NA, max_len - length(i)))))

#lengths of all sequences
seq_lengths <- data.frame(len = c(lengths(seq_pos), lengths(seq_neg)),
                          tar = c(rep("pos", length(seq_pos)),
                                 rep("neg", length(seq_neg))))

save(seq_lengths, file = "gcb_abstract_poster/poster_data.RData")


library(ggplot2)
ggplot(filter(seq_lengths, len < 50), aes(x = len, fill = tar)) +
  geom_density(alpha = 0.5)
