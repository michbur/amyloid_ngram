library(seqinr)
library(randomForest)
library(dplyr)
library(biogram)

source("aa_encodings2.R") #creates encodings for data

#load comparison of rep33 and AmyLoad to remove doubles
compr <- read.csv2("rep33vsAmyLoad.csv", row.names = NULL, stringsAsFactors = FALSE)

#sequences of amyloids existing both in AmyLoad and rep33
doubles <- filter(compr, AmyLoad != "brak") %>% 
  #mutate(nice_name = sapply(strsplit(name, "_", fixed = TRUE), first)) %>%
  select(AmyLoad) %>% unlist(use.names = FALSE) %>% unique

seq_pos <- read.fasta("gcb_abstract_poster/amyloid_pos_full.fasta", seqtype = "AA")

pos_pasted <- sapply(seq_pos, function(i) paste0(i, collapse = ""))

#which sequences should be removed from the positive data set
pos_remove <- which(pos_pasted %in% doubles)

pos_train <- seq_pos[-pos_remove]
neg_train <- read.fasta("gcb_abstract_poster/amyloid_neg_full.fasta", seqtype = "AA")

filter_length <- function(seq, max_length) {
  seq[lengths(seq) <= max_length & lengths(seq) > 5]
}

create_gl <- function(seq)
  lapply(1L:nrow(seq), function(i) {
    res <- do.call(rbind, strsplit(decode_ngrams(seq2ngrams(seq[i, ][!is.na(seq[i, ])], 6, a()[-1])), ""))
    cbind(res, id = paste0("P", rep(i, nrow(res))))
  })

get_bitrigrams <- function(seq, aa_group) 
  lapply(seq, function(single_protein) {
    bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                             ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                             seq = degenerate(single_protein[, -7], aa_group),
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

pos_classifier <- function(pos_dat, neg_dat, aa_group) {
  pos_fil <- filter_length(pos_dat, 6)
  neg_fil <- filter_length(neg_dat, 6)
  
  all_ft <- do.call(rbind, flist2matrix(c(pos_fil, neg_fil)) %>% 
                      create_gl() %>% get_bitrigrams(aa_group = aa_group))
  
  test_bis <- test_features(target = c(rep(1, length(pos_fil)), rep(0, length(neg_fil))),
                            features = all_ft, adjust = NULL)
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  
  list(model = randomForest(x = all_ft[, imp_bigrams], as.factor(c(rep("pos", length(pos_fil)), 
                                                                   rep("neg", length(neg_fil))))),
       imps <- imp_bigrams)
}

neg_classifier <- function(pos_dat, neg_dat, aa_group) {
  pos_fil <- filter_length(pos_dat, 15)
  neg_fil <- filter_length(neg_dat, 15)
  
  all_ft <- do.call(rbind, flist2matrix(c(pos_fil, neg_fil)) %>% 
                      create_gl() %>% get_bitrigrams(aa_group = aa_group))
  
  test_bis <- test_features(target = c(rep(1, length(pos_fil)), rep(0, length(neg_fil))),
                            features = all_ft, adjust = NULL)
  imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
  
  list(model = randomForest(x = all_ft[, imp_bigrams], as.factor(c(rep("pos", length(pos_fil)), 
                                                                   rep("neg", length(neg_fil))))),
       imps <- imp_bigrams)
}

pos_rf <- pos_classifier(pos_train, neg_train, aa_groups[[45]])
neg_rf <- neg_classifier(pos_train, neg_train, aa_groups[[87]])

save(pos_rf, neg_rf, file = "tmp_res.RData")