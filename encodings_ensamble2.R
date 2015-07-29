#training: we extract 5-grams from sequences. All 5-grams are labelled as the sequences of origin. 
#In case of the 15-aa-long amyloids, we may introduce a lot of false positive n-grams
#moreover, training only on 6-mers, we have only aroung 50% of all non-amyloid sequences. The negative data set may be too small
#hypothesis:
#The factors above may make classifier trained on 6-mer weak in recognising negative cases and trained on 15-mer 
#weak in recognising positive cases

source("read_data.R")


library(cvTools)
library(randomForest)
library(pbapply)
library(biogram)

load("encoded_amyloids2.RData")

pos_constant <- apply(pos_data, 1, function(i) length(i) - sum(is.na(i))) <= 6
neg_constant <- apply(neg_data, 1, function(i) length(i) - sum(is.na(i))) > 15

fold_list <- lapply(list(pos_constant, !pos_constant, neg_constant, !neg_constant), function(single_n) {
  folded <- cvFolds(sum(single_n), K = 5)
  cbind(id = which(single_n)[folded[["subsets"]]], which = folded[["which"]])
})

res <- lapply(1L:5, function(fold) {
  lapply(c(86, 91), function(group_id) {
    train_pos <- do.call(rbind, coded_pos[[group_id]][fold_list[[1]][fold_list[[1]][, "which"] != fold, "id"]])
    train_neg <- do.call(rbind, coded_neg[[group_id]][fold_list[[3]][fold_list[[3]][, "which"] != fold, "id"]])
    
    test_pos <- rbind(do.call(rbind, coded_pos[[group_id]][fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"]]),
                      do.call(rbind, coded_pos[[group_id]][fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"]]))
    test_neg <- rbind(do.call(rbind, coded_neg[[group_id]][fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"]]), 
                      do.call(rbind, coded_neg[[group_id]][fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"]]))
    
    test_bis <- test_features(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg))),
                              rbind(train_pos, train_neg), adjust = NULL)
    imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
    
    
    model_cv <- randomForest(x = rbind(train_pos, train_neg)[, imp_bigrams], 
                             as.factor(c(rep(1, nrow(train_pos)), rep(0, nrow(train_neg)))))
    
    #number of n-grams from protein
    ngram_prots_pos <- c(sapply(coded_pos[[group_id]][fold_list[[1]][fold_list[[1]][, "which"] == fold, "id"]], nrow),
                         sapply(coded_pos[[group_id]][fold_list[[2]][fold_list[[2]][, "which"] == fold, "id"]], nrow))
    ngram_prots_neg <- c(sapply(coded_neg[[group_id]][fold_list[[3]][fold_list[[3]][, "which"] == fold, "id"]], nrow),
                         sapply(coded_neg[[group_id]][fold_list[[4]][fold_list[[4]][, "which"] == fold, "id"]], nrow))
    
    preds <- cbind(predict(model_cv, rbind(test_pos, test_neg), type = "prob")[, 2], 
                   c(rep(1, nrow(test_pos)), rep(0, nrow(test_neg))),
                   c(unlist(lapply(1L:length(ngram_prots_pos), function(prot_id)
                     rep(prot_id, ngram_prots_pos[prot_id]))), 
                     unlist(lapply(1L:length(ngram_prots_neg), function(prot_id)
                       rep(prot_id, ngram_prots_neg[prot_id])))))
    
    
    list(preds, imp_bigrams)
  })
})

library(hmeasure)
library(dplyr)

fold_cv <- lapply(res, function(single_fold) {
  enc_cv <- sapply(single_fold, function(single_enc) {
    dat <- data.frame(single_enc[[1]])
    dat1 <- dat %>% filter(X2 == 1) %>% group_by(X3) %>% summarize(pred = max(X1))
    dat0 <- dat %>% filter(X2 == 0) %>% group_by(X3) %>% summarize(pred = max(X1))
    real_labels <- c(rep(1, nrow(dat1)), rep(0, nrow(dat0)))
    unlist(HMeasure(real_labels, c(dat1[["pred"]], dat0[["pred"]]))[["metrics"]])
  })
  cbind(enc = 1L:ncol(enc_cv), t(enc_cv))
})

fold_cv <- do.call(rbind, lapply(1L:length(fold_cv), function(single_fold_id)
  cbind(fold = rep(single_fold_id, nrow(fold_cv[[single_fold_id]])), fold_cv[[single_fold_id]])))

data.frame(fold_cv) %>% group_by(enc) %>% summarize(mean(AUC), mean(Sens), mean(Spec))
#enc mean(AUC) mean(Sens) mean(Spec)
#1   1 0.8275416  0.7448421  0.7655226
#2   2 0.8334615  0.8166667  0.7234322