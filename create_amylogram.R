library(seqinr)
library(randomForest)
library(dplyr)
library(biogram)
library(hmeasure)

source("aa_encodings2.R") #creates encodings for data

seq_pos <- read.fasta("gcb_abstract_poster/amyloid_pos_full.fasta", seqtype = "AA")
seq_neg <- read.fasta("gcb_abstract_poster/amyloid_neg_full.fasta", seqtype = "AA")

filter_length <- function(seq, max_length) {
  seq[lengths(seq) <= max_length & lengths(seq) > 5]
}

create_gl <- function(seq)
  lapply(1L:nrow(seq), function(i) {
    res <- do.call(rbind, strsplit(decode_ngrams(seq2ngrams(seq[i, ][!is.na(seq[i, ])], min_subseq_length, a()[-1])), ""))
    cbind(res, id = paste0("P", rep(i, nrow(res))))
  })

get_bitrigrams <- function(seq, aa_group) 
  lapply(seq, function(single_protein) {
    bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                             ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                             seq = degenerate(single_protein[, -(min_subseq_length + 1)], aa_group),
                                             u = as.character(1L:length(aa_group))))
    
    bitrigrams <- bitrigrams > 0
    storage.mode(bitrigrams) <- "integer"
    
    bitrigrams
  })


# fasta list to matrix
flist2matrix <- function(x) {
  max_len <- max(lengths(x))
  
  t(sapply(x, function(i)
    c(i, rep(NA, max_len - length(i)))))
}

make_classifier <- function(pos_dat, neg_dat, aa_group, max_lenth) {
  pos_fil <- filter_length(pos_dat, max_lenth)
  neg_fil <- filter_length(neg_dat, max_lenth)
  
  pos_ft <- do.call(rbind, flist2matrix(pos_fil) %>% 
                      create_gl() %>% 
                      get_bitrigrams(aa_group = aa_group))
  
  neg_ft <- do.call(rbind, flist2matrix(neg_fil) %>% 
                      create_gl() %>% 
                      get_bitrigrams(aa_group = aa_group))
  all_ft <- rbind(pos_ft, neg_ft)
  
  test_bis <- test_features(target = c(rep(1, nrow(pos_ft)), rep(0, nrow(neg_ft))),
                            features = all_ft, adjust = NULL)
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  
  list(model = randomForest(x = all_ft[, imp_bigrams], as.factor(c(rep("pos", nrow(pos_ft)), 
                                                                   rep("neg", nrow(neg_ft))))),
       imps = imp_bigrams)
}



pos_rf <- make_classifier(seq_pos, seq_neg, aa_groups[[45]], 6)
neg_rf <- make_classifier(seq_pos, seq_neg, aa_groups[[87]], 15)