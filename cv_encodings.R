source("read_data.R")

source("aa_encodings.R") #creates encodings for data

source("encoded_amyloids.R") #codes amyloids and saves them

library(cvTools)
library(randomForest)
library(pbapply)


fold_res <- pblapply(1L:4, function(dummy) {
  pos_ids <- cvFolds(sum(targets == 1), K = 5)
  neg_ids <- cvFolds(sum(targets != 1), K = 5)
  lapply(1L:5, function(fold) {
    lapply(1L:length(aa_groups), function(group_id) {
      train_pos <- as.matrix(coded_pos[[group_id]])[pos_ids[[5]] != fold, ]
      train_neg <- as.matrix(coded_neg[[group_id]])[neg_ids[[5]] != fold, ]
      
      test_pos <- as.matrix(coded_pos[[group_id]])[pos_ids[[5]] == fold, ]
      test_neg <- as.matrix(coded_neg[[group_id]])[neg_ids[[5]] == fold, ]
      
      test_bis <- test_features(c(rep(1, sum(pos_ids[[5]] != fold)), rep(0, sum(neg_ids[[5]] != fold))),
                                rbind(train_pos, train_neg), adjust = NULL)
      imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
      
      
      model_cv <- randomForest(x = rbind(train_pos, train_neg)[, imp_bigrams], 
                               as.factor(c(rep(1, sum(pos_ids[[5]] != fold)), rep(0, sum(neg_ids[[5]] != fold)))))
      
      preds <- cbind(as.numeric(as.character(predict(model_cv, rbind(test_pos, test_neg)))), 
                     c(rep(1, sum(pos_ids[[5]] == fold)), rep(0, sum(neg_ids[[5]] == fold))))
      preds
    })
  })
})

fold_res6 <- pblapply(1L:4, function(dummy) {
  
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
                                rbind(train_pos, train_neg), adjust = NULL)
      imp_bigrams <- cut(test_bis, breaks = c(0, 0.1, 1))[[1]]
      
      
      model_cv <- randomForest(x = rbind(train_pos, train_neg)[, imp_bigrams], 
                               as.factor(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg)))))
      
      preds <- cbind(predict(model_cv, rbind(test_pos, test_neg), type = "prob")[, 2], 
                     c(rep(1, nrow(test_pos)), rep(0, nrow(test_neg))))
      preds
    })
  })
})

save(fold_res, fold_res6, file = "amyloid_fold_res.RData")
