source("read_data.R")

source("aa_encodings.R")

library(biogram)
library(cvTools)
library(randomForest)
library(pbapply)

pos_data <- all_seqs[as.logical(targets), ]
neg_data <- all_seqs[!as.logical(targets), ]

seq_lengths <- apply(pos_data, 1, function(i) length(i) - sum(is.na(i)))
pos_6 <- sapply(seq_pos, length) == 6
neg_6 <- sapply(seq_neg, length) != 6


coded_pos <- lapply(aa_groups, function(single_group) {
  bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                           ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                           seq = degenerate(as.matrix(pos_data), single_group),
                                           u = as.character(1L:length(single_group))))
  
  bitrigrams <- bitrigrams > 0
  storage.mode(bitrigrams) <- "integer"
  
  bitrigrams
})


coded_neg <- lapply(aa_groups, function(single_group) {
  bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                           ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                           seq = degenerate(as.matrix(neg_data), single_group),
                                           u = as.character(1L:length(single_group))))
  
  bitrigrams <- bitrigrams > 0
  storage.mode(bitrigrams) <- "integer"
  
  bitrigrams
})


fold_res6 <- pblapply(1L:50, function(dummy) {
  
  fold_list <- lapply(list(pos_6, !pos_6, neg_6, !neg_6), function(single_n) {
    folded <- cvFolds(sum(single_n), K = 5)
    cbind(id = which(single_n)[folded[["subsets"]]], which = folded[["which"]])
  })
  
  lapply(1L:5, function(fold) {
    lapply(1L:length(aa_groups), function(group_id) {
      train_pos <- as.matrix(coded_pos[[group_id]])[fold_list[[1]][fold_list[[1]][, "which"] != fold, "id"], ]
      train_neg <- as.matrix(coded_neg[[group_id]])[fold_list[[3]][fold_list[[3]][, "which"] != fold, "id"], ]
      
      test_pos <- rbind(as.matrix(coded_pos[[group_id]])[fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"], ],
                        as.matrix(coded_pos[[group_id]])[fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"], ])
      test_neg <- rbind(as.matrix(coded_neg[[group_id]])[fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"], ], 
                        as.matrix(coded_neg[[group_id]])[fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"], ])
      
      test_bis <- test_features(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg))),
                                rbind(train_pos, train_neg))
      imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
      
      
      model_cv <- randomForest(x = rbind(train_pos, train_neg)[, imp_bigrams], 
                               as.factor(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg)))))
      
      preds <- cbind(predict(model_cv, rbind(test_pos, test_neg), type = "prob")[, 2], 
                     c(rep(1, nrow(test_pos)), rep(0, nrow(test_neg))))
      preds
    })
  })
})


save(fold_res6, file = "amyloid_fold_res.RData")
