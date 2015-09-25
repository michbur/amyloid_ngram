library(seqinr)
library(dplyr)
library(biogram)
library(reshape2)
library(ggplot2)

source("aa_encodings2.R") #creates encodings for data

spec = aa_groups[[45]]
sens = aa_groups[[87]]

seq_pos <- read.fasta("gcb_abstract_poster/amyloid_pos_full.fasta", seqtype = "AA")
seq_neg <- read.fasta("gcb_abstract_poster/amyloid_neg_full.fasta", seqtype = "AA")

flist2matrix <- function(x) {
  max_len <- max(lengths(x))
  
  t(sapply(x, function(i)
    c(i, rep(NA, max_len - length(i)))))
}

get_freqs <- function(aa_group_id, n = 1) {
  rbind(flist2matrix(seq_pos) %>% 
                    degenerate(element_groups = aa_group_id) %>%
                    count_ngrams(n = n, u = 1L:length(aa_group_id)) %>%
                    as.matrix %>%
                    data.frame %>%
                    cbind(., tar = "yes"),
                  flist2matrix(seq_neg) %>% 
                    degenerate(element_groups = aa_group_id) %>%
                    count_ngrams(n = n, u = 1L:length(aa_group_id)) %>%
                    as.matrix %>%
                    data.frame %>%
                    cbind(., tar = "no")) %>%
    melt %>% group_by(tar, variable) %>% summarise(freq = sum(value)) %>% ungroup %>% 
    group_by(tar) %>% mutate(freq = freq/sum(freq)) %>% data.frame
}

amylo_freq <- rbind(cbind(get_freqs(aa_groups[[45]]), enc = "Best specificity", n = 1),
                    cbind(get_freqs(aa_groups[[87]]), enc = "Best sensitivity", n = 1),
                    cbind(get_freqs(aa_groups[[45]], n= 2), enc = "Best specificity", n = 2),
                    cbind(get_freqs(aa_groups[[87]], n= 2), enc = "Best sensitivity", n = 2))

levels(amylo_freq[["variable"]]) <- decode_ngrams(sapply(levels(amylo_freq[["variable"]]), function(i) substr(i, 2, nchar(i))))


save(amylo_freq, spec, sens, file = "./gcb_abstract_poster/specsens.RData")

ggplot(amylo_freq, aes(x = variable, y = freq, fill = tar, colour = tar)) +
  geom_bar(stat = "identity", position = "dodge") + 
  scale_y_continuous("Frequency") +
  scale_x_discrete("Group ID\n") +
  scale_fill_discrete("Amyloid") +
  facet_wrap( ~ enc, scales = "free") +
  guides(colour = FALSE) 