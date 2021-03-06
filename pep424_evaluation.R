library(seqinr)
library(randomForest)
library(dplyr)
library(biogram)
library(hmeasure)

source("aa_encodings2.R") #creates encodings for data

pep424_tab <- read.table("pep424_evaluation.txt", sep = "\t")
levels(pep424_tab[[3]]) <- c(0, 0, 1, 1)


pep424_prots <- lapply(as.character(pep424_tab[[2]]) %>% strsplit(split = ""), function(i)
  i[i != " "])

# write pasta files with seqs
write.fasta(pep424_prots[lengths(pep424_prots) > 5], 
            names = as.character(pep424_tab[[1]]), file.out = "pep424.fasta")

write.fasta(pep424_prots[lengths(pep424_prots) > 5], 
            names = paste0("prot", 1L:sum(lengths(pep424_prots) > 5)), file.out = "pep424_better_names.fasta")

#reads data and removes seqs already present in test set
get_train_seqs <- function(file_name, pattern_vector) {
  seq <- read.fasta(file_name, seqtype = "AA")
  
  seq_pasted <- sapply(seq, function(i) paste0(i, collapse = ""))
  
  seq_duplicated <- lapply(seq_pasted, function(i)
    grep(i, pattern_vector)) %>% unlist %>% unique
  
  seq[-seq_duplicated]
}

seq_pos <- get_train_seqs("gcb_abstract_poster/amyloid_pos_full.fasta", sapply(pep424_prots, paste0, collapse = ""))
seq_neg <- get_train_seqs("gcb_abstract_poster/amyloid_neg_full.fasta", sapply(pep424_prots, paste0, collapse = ""))

min_subseq_length <- 5

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

test_pos_data <- pep424_prots[lengths(pep424_prots) > min_subseq_length] %>%
  flist2matrix %>%
  create_gl %>% 
  get_bitrigrams(aa_group = aa_groups[[45]]) %>%
  do.call(rbind, .) 

pos_pred <- predict(pos_rf[[1]], test_pos_data[, pos_rf[[2]]], type = "prob")[, 2]

test_neg_data <- pep424_prots[lengths(pep424_prots) > min_subseq_length] %>%
  flist2matrix %>%
  create_gl %>% 
  get_bitrigrams(aa_group = aa_groups[[87]]) %>%
  do.call(rbind, .) 

neg_pred <- predict(neg_rf[[1]], test_neg_data[, neg_rf[[2]]], type = "prob")[, 2]

final_pred <- (pos_pred + neg_pred)/2

long_lengths <- lengths(pep424_prots)[lengths(pep424_prots) > min_subseq_length]
prot_preds <- data.frame(id = unlist(lapply(1L:length(long_lengths), function(i) 
  rep(i, long_lengths[i] - min_subseq_length + 1))), pred = final_pred) %>% 
  group_by(id) %>% summarise(pred = max(pred)) %>% select(pred) %>%
  unlist

#real labels
real_labels <- as.numeric(as.character(pep424_tab[[3]]))[lengths(pep424_prots) > min_subseq_length]


#fold amyloid preds
foldAmyloid_data <- readLines("FoldAmyloid_pred.txt")
foldAmyloid_preds <- grepl("1", 
                           foldAmyloid_data[grep("Options: Scale = Expected number of contacts 8A", 
                                                 foldAmyloid_data) + 1])

#PASTA 2.0 preds
file_names <- list.files("./pasta2_preds/")[grepl("fasta.seq.best_pairings_list.dat", list.files("./pasta2_preds/"))]

pasta2_prot_id <- as.numeric(sapply(strsplit(file_names, ".", fixed = TRUE), function(i) substr(i[1], 5, nchar(i[1]))))

pasta2_preds <- sapply(file_names, function(i) {
  all_lines <- readLines(paste0("./pasta2_preds/", i))
  as.numeric(strsplit(strsplit(all_lines, "PASTA energy ")[[1]][2], "*[ ]")[[1]][1])
}) < -4



metrics <- HMeasure(real_labels, data.frame(AmyloGram = prot_preds,
                                            PASTA2 = pasta2_preds[order(pasta2_prot_id)],
                                            FoldAmyloid = foldAmyloid_preds))[["metrics"]][, c("AUC", "Sens", "Spec")]

save(metrics, file = "gcb_abstract_poster/metrics.RData")
