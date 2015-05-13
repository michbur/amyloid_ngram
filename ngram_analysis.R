#goals:
#significant 2- and 3-grams (all, hexamers, other than hexamers)
#the length of the signal

source("read_data.R")
library(biogram)
library(dplyr)

#lengths of all sequences
seq_lengths <- ncol(all_seqs) - apply(all_seqs, 1, function(i) sum(is.na(i)))

#all data
distances1 <- 0
distances2 <- 0L:3
distances3 <- unlist(apply(expand.grid(0L:1, 0L:1), 1, list), recursive = FALSE)
distances <- c(list(distances1), as.list(distances2), distances3)

sizes <- c(1, rep(2, length(distances2)), rep(3, length(distances3)))

bitrigrams <- count_multigrams(seq = all_seqs, 
                               ns = sizes, 
                               u = tolower(a()[-1]), 
                               ds = distances,
                               pos = FALSE)
bitrigrams_bin <- as.matrix(bitrigrams) > 0
storage.mode(bitrigrams_bin) <- "integer"
test_res <- test_features(targets, bitrigrams_bin, adjust = NULL)
imp_ngrams <- na.omit(cut(test_res, breaks = c(0, 1e-3, 1))[[1]])

n_l <- sapply(strsplit(ngrams2df(imp_ngrams)[, "ngram"], ".", fixed = TRUE), length)

ngram_tables <- lapply(unique(n_l), function(i) table_ngrams(all_seqs, imp_ngrams[n_l == i], targets))


total_ngram <- sapply(1L:length(sizes), function(id_size)
  sapply(sort(unique(seq_lengths)), function(peptide_length)
    length(get_ngrams_ind(peptide_length, sizes[id_size], distances[[id_size]])[[1]])))

colnames(total_ngram) <- sizes
rownames(total_ngram) <- sort(unique(seq_lengths))

#total number of positions in positve n-grams
total_ngrams_pos <- rowSums(sapply(seq_lengths[targets == 1], function(i)
  total_ngram[as.character(i), ]))

#total number of positions in negative n-grams
total_ngrams_neg <- rowSums(sapply(seq_lengths[targets == 0], function(i)
  total_ngram[as.character(i), ]))

ngram_freqs <- do.call(rbind, ngram_tables)

#size of important n-grams
imp_size <- ngram_freqs[, "ngram"] %>% 
  as.character %>% 
  ngrams2df %>% 
  select(ngram) %>%
  unlist %>%
  strsplit(".", fixed = TRUE) %>%
  lapply(length) %>%
  unlist

#distance of important n-grams
imp_dist <- ngram_freqs[, "ngram"] %>% 
  as.character %>% 
  ngrams2df %>% 
  select(distance) %>%
  unlist 



pasted_dists <- lapply(distances, function(i) paste(i, collapse = "."))

all_sizes_anal <- cbind(ngram_freqs, 
      do.call(rbind, lapply(1L:length(sizes), function(i)
        data.frame(ftarget0 = ngram_freqs[imp_size == sizes[i] & imp_dist == pasted_dists[i], "target0"]/total_ngrams_neg[i],
                   ftarget1 = ngram_freqs[imp_size == sizes[i] & imp_dist == pasted_dists[i], "target1"]/total_ngrams_pos[i]
        ))),
      decoded = ngram_freqs[["ngram"]] %>% as.character %>% decode_ngrams
)

bigger6_anal <- analyze_ngrams(all_seqs[seq_lengths > 6, ], targets[seq_lengths > 6])
only6_anal <- analyze_ngrams(all_seqs[seq_lengths == 6, ], targets[seq_lengths == 6])

write.csv2(all_sizes_anal, "all_sizes.csv")
write.csv2(bigger6_anal, "bigger6.csv")
write.csv2(only6_anal, "only6.csv")
