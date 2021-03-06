library(seqinr)
library(randomForest)
library(dplyr)
library(biogram)
library(hmeasure)

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

#hexamers here
min_subseq_length <- 6

filter_length <- function(seq, max_length) {
  seq[lengths(seq) <= max_length & lengths(seq) > min_subseq_length - 1]
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



# read test data --------------------------------

all_reg <- readLines("reg33_full.txt")
all_reg_id <- matrix(1L:length(all_reg), ncol = 3, byrow = TRUE)


real_labels <- do.call(rbind, lapply(1L:nrow(all_reg_id), function(i) {
  prot_seq <- all_reg[all_reg_id[i, 2]] %>%
    strsplit(split = "") %>%
    unlist
  amyl <- rep(0, length(prot_seq))
  
  amyl_regs <- all_reg[all_reg_id[i, 3]] %>%
    strsplit(split = " ") %>%
    unlist %>% 
    as.numeric %>%
    matrix(nrow = 2) %>%
    apply(2, function(j) j[1]:j[2] + 1) %>%
    unlist %>%
    as.vector
  #-1 because reg33 counts 0, 1
  amyl[amyl_regs - 1] <- 1
  
  res <- data.frame(id = rep(i, length(prot_seq)), prot_seq = prot_seq, amyl = amyl)
  colnames(res) <- c("id", "prot_seq", "amyl")
  res
}))


test_pos_data <- all_reg[all_reg_id[, 2]] %>%
  strsplit(split = "") %>%
  flist2matrix %>%
  create_gl %>% 
  get_bitrigrams(aa_group = aa_groups[[45]]) %>%
  do.call(rbind, .) 

pos_pred <- predict(pos_rf[[1]], test_pos_data[, pos_rf[[2]]], type = "prob")[, 2]

test_neg_data <- all_reg[all_reg_id[, 2]] %>%
  strsplit(split = "") %>%
  flist2matrix %>%
  create_gl %>% 
  get_bitrigrams(aa_group = aa_groups[[87]]) %>%
  do.call(rbind, .) 

neg_pred <- predict(neg_rf[[1]], test_neg_data[, neg_rf[[2]]], type = "prob")[, 2]

mean_pred <- all_reg[all_reg_id[, 2]] %>%
  strsplit(split = "") %>%
  flist2matrix %>%
  create_gl %>%
  do.call(rbind, .) %>%
  data.frame %>%
  select(id) %>%
  cbind(., pred = (neg_pred + pos_pred)/2)



biogram_res <- unlist(lapply(levels(mean_pred[["id"]]), function(single_protein) {
  hexapreds <- filter(mean_pred, id == single_protein) %>% 
    select(pred) %>% 
    unlist 
  
  sapply(1L:(length(hexapreds) + 5), function(i) {
    range_low <- ifelse(i - 5 < 1, 1, i - 5)
    range_up <- ifelse(i + 5 > length(hexapreds), length(hexapreds), i + 5)
    print(max(hexapreds[range_low:range_up]))
    max(hexapreds[range_low:range_up])
  })
}))

bench_res <- cbind(real_labels, pred = biogram_res)

HMeasure(bench_res[["amyl"]], bench_res[["pred"]])[["metrics"]]
