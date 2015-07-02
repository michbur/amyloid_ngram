#goals:
#significant 2- and 3-grams (all, hexamers, other than hexamers)
#the length of the signal

source("read_data.R")
source("analyze_ngrams_function.R")

library(biogram)
library(dplyr)

all_sizes_anal <- analyze_ngrams(all_seqs, targets)
bigger6_anal <- analyze_ngrams(all_seqs[seq_lengths > 6, ], targets[seq_lengths > 6])
only6_anal <- analyze_ngrams(all_seqs[seq_lengths == 6, ], targets[seq_lengths == 6])

write.csv2(all_sizes_anal, "all_sizes.csv")
write.csv2(bigger6_anal, "bigger6.csv")
write.csv2(only6_anal, "only6.csv")
