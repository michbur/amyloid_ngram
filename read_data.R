library(seqinr)

seq_pos <- read.fasta("amyloid_pos.fasta", seqtype = "AA")
seq_neg <- read.fasta("amyloid_neg.fasta", seqtype = "AA")

max_len <- max(sapply(c(seq_pos, seq_neg), length))

all_seqs <- t(sapply(c(seq_pos, seq_neg), function(i)
  c(i, rep(NA, max_len - length(i)))))

#lengths of all sequences
seq_lengths <- ncol(all_seqs) - apply(all_seqs, 1, function(i) sum(is.na(i)))

#remove too short sequences
all_seqs <- all_seqs[seq_lengths > 4, ]

targets <- c(rep(1, length(seq_pos) -1), rep(0, length(seq_neg)))
