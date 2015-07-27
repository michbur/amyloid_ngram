library(seqinr)
library(dplyr)

seq_pos <- read.fasta("gcb_abstract_poster/amyloid_pos_full.fasta", seqtype = "AA")
seq_neg <- read.fasta("gcb_abstract_poster/amyloid_neg_full.fasta", seqtype = "AA")

max_len <- max(sapply(c(seq_pos, seq_neg), length))

all_seqs <- t(sapply(c(seq_pos, seq_neg), function(i)
  c(i, rep(NA, max_len - length(i)))))

#lengths of all sequences
seq_lengths <- data.frame(len = c(lengths(seq_pos), lengths(seq_neg)),
                          tar = c(rep("yes", length(seq_pos)),
                                 rep("no", length(seq_neg))))
#sequence, length in the invervals
seq_inter <- data.frame(len = cut(seq_lengths[["len"]], 
                                  c(min(seq_lengths[["len"]]), 5, 6, 10, 15, max(seq_lengths[["len"]])), 
                                  include.lowest = TRUE),
                        tar = seq_lengths[["tar"]])

#sequence, length in the invervals, crosstable
seq_inter_ct <- data.frame(xtabs(~ len + tar, seq_inter)) %>% group_by(tar) %>% mutate(prop = Freq/sum(Freq)) %>% ungroup

save(seq_inter_ct, file = "gcb_abstract_poster/poster_data.RData")

# Density plot that may be useful
# seq_lengthsf <- filter(seq_lengths, len < 41)
# seq_lengthsf[["len"]] <- factor(seq_lengthsf[["len"]], levels = 4L:40)
# 
# ggplot(data.frame(table(seq_lengthsf[["len"]], seq_lengthsf[["tar"]])), aes(x = Var1, y = Freq, fill = Var2)) +
#   geom_bar(alpha = 0.5, position = "dodge", stat = "identity") +
#   scale_fill_discrete("Amyloid") +
#   scale_x_discrete("Length") +
#   scale_y_continuous("Count")
