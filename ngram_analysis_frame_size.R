library(hmeasure)
library(reshape2)
library(seqinr)
#library(ggvis)
library(dplyr)
library(magrittr)
library(ggplot2)
library(gridExtra)

source("aa_encodings2.R")
load("amyloid_fold_res2.RData")

names(fold_res[[1]])
head(fold_res[[1]][[1]][[1]][[1]][[1]])

len_cv <- lapply(fold_res, function(single_length) {
  rep_cv <- lapply(single_length, function(single_rep) {
    fold_cv <- lapply(single_rep, function(single_fold) {
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
  })
  library(hmeasure)
  library(reshape2)
  library(seqinr)
  #library(ggvis)
  library(dplyr)
  library(magrittr)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  
  source("aa_encodings2.R")
  load("amyloid_fold_res2.RData")
  
  my_max <- function(x){
    w <- NULL
    if(length(x)>1){
      w <- mean(sort(x,TRUE)[1:2])
    } else {
      w <- x
    }
    w
  }
  
  max2 <- function(x){
    if(length(x)>1){
      return(sort(x,TRUE)[2])
    } else {
      return(x)
    }
  }
  
  calc_mcc <- function(metrics) {
    TP <- as.numeric(metrics[["TP"]])
    TN <- metrics[["TN"]]
    FP <- metrics[["FP"]]
    FN <- metrics[["FN"]]
    data.frame(metrics, mcc = unname((TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))))
  }
  
  len_cv <- lapply(fold_res, function(single_length) {
    rep_cv <- lapply(single_length, function(single_rep) {
      fold_cv <- lapply(single_rep, function(single_fold) {
        enc_cv <- sapply(single_fold, function(single_enc) {
          dat <- data.frame(single_enc[[1]])
          dat1 <- dat %>% filter(X2 == 1) %>% group_by(X3) %>% summarize(pred = max(X1))
          dat0 <- dat %>% filter(X2 == 0) %>% group_by(X3) %>% summarize(pred = max(X1))
          real_labels <- c(rep(1, nrow(dat1)), rep(0, nrow(dat0)))
          unlist(calc_mcc(HMeasure(real_labels, c(dat1[["pred"]], dat0[["pred"]]))[["metrics"]]))
        })
        cbind(enc = 1L:ncol(enc_cv), t(enc_cv))
      })
      
      fold_cv <- do.call(rbind, lapply(1L:length(fold_cv), function(single_fold_id)
        cbind(fold = rep(single_fold_id, nrow(fold_cv[[single_fold_id]])), fold_cv[[single_fold_id]])))
    })
    
    rep_cv <- do.call(rbind, lapply(1L:length(rep_cv), function(single_rep_id)
      cbind(repetition = rep(single_rep_id, nrow(rep_cv[[single_rep_id]])), rep_cv[[single_rep_id]])))
  })
  len_cv2 <- do.call(rbind, lapply(1L:length(len_cv), function(single_len_id)
    cbind(len = rep(single_len_id, nrow(len_cv[[single_len_id]])), len_cv[[single_len_id]])))
  
  len_cv2 <- data.frame(len_cv2)
  len_cv2[["len"]] <- as.factor(len_cv2[["len"]])
  levels(len_cv2[["len"]]) <- c("6", "<11", "<16")
  len_cv2[["repetition"]] <- as.factor(len_cv2[["repetition"]])
  len_cv2[["enc"]] <- as.factor(len_cv2[["enc"]])
  len_cv2[["fold"]] <- as.factor(len_cv2[["fold"]])
  
  
  res <- len_cv2 %>% group_by(len, enc, repetition, fold) %>% summarize(mAUC = mean(AUC), mSens = mean(Sens), 
                                                                        mSpec = mean(Spec), mMcc = mean(mcc)) %>%
    group_by(len, enc) %>% summarize(mAUC = mean(mAUC), mSens = mean(mSens), mSpec = mean(mSpec), mMcc = mean(mMcc)) %>%
    filter(mAUC > quantile(mAUC, 0.95)) 
  res <- res[order(res[["mAUC"]], decreasing = TRUE), ]
  
  len_cv <- lapply(fold_res, function(single_length) {
    rep_cv <- lapply(single_length, function(single_rep) {
      fold_cv <- lapply(single_rep, function(single_fold) {
        enc_cv <- sapply(single_fold, function(single_enc) {
          dat <- data.frame(single_enc[[1]])
          dat1 <- dat %>% filter(X2 == 1) %>% group_by(X3) %>% summarize(pred = my_max(X1))
          dat0 <- dat %>% filter(X2 == 0) %>% group_by(X3) %>% summarize(pred = my_max(X1))
          real_labels <- c(rep(1, nrow(dat1)), rep(0, nrow(dat0)))
          unlist(calc_mcc(HMeasure(real_labels, c(dat1[["pred"]], dat0[["pred"]]))[["metrics"]]))
        })
        cbind(enc = 1L:ncol(enc_cv), t(enc_cv))
      })
      
      fold_cv <- do.call(rbind, lapply(1L:length(fold_cv), function(single_fold_id)
        cbind(fold = rep(single_fold_id, nrow(fold_cv[[single_fold_id]])), fold_cv[[single_fold_id]])))
    })
    
    rep_cv <- do.call(rbind, lapply(1L:length(rep_cv), function(single_rep_id)
      cbind(repetition = rep(single_rep_id, nrow(rep_cv[[single_rep_id]])), rep_cv[[single_rep_id]])))
  })
  len_cv2 <- do.call(rbind, lapply(1L:length(len_cv), function(single_len_id)
    cbind(len = rep(single_len_id, nrow(len_cv[[single_len_id]])), len_cv[[single_len_id]])))
  
  len_cv2 <- data.frame(len_cv2)
  len_cv2[["len"]] <- as.factor(len_cv2[["len"]])
  levels(len_cv2[["len"]]) <- c("6", "<11", "<16")
  len_cv2[["repetition"]] <- as.factor(len_cv2[["repetition"]])
  len_cv2[["enc"]] <- as.factor(len_cv2[["enc"]])
  len_cv2[["fold"]] <- as.factor(len_cv2[["fold"]])
  
  
  res <- len_cv2 %>% group_by(len, enc, repetition, fold) %>% summarize(mAUC = mean(AUC), mSens = mean(Sens), 
                                                                        mSpec = mean(Spec), mMcc = mean(mcc)) %>%
    group_by(len, enc) %>% summarize(mAUC = mean(mAUC), mSens = mean(mSens), mSpec = mean(mSpec), mMcc = mean(mMcc)) %>%
    filter(mAUC > quantile(mAUC, 0.95)) 
  res2 <- res[order(res[["mAUC"]], decreasing = TRUE), ]
  
  len_cv <- lapply(fold_res, function(single_length) {
    rep_cv <- lapply(single_length, function(single_rep) {
      fold_cv <- lapply(single_rep, function(single_fold) {
        enc_cv <- sapply(single_fold, function(single_enc) {
          dat <- data.frame(single_enc[[1]])
          dat1 <- dat %>% filter(X2 == 1) %>% group_by(X3) %>% summarize(pred = max2(X1))
          dat0 <- dat %>% filter(X2 == 0) %>% group_by(X3) %>% summarize(pred = max2(X1))
          real_labels <- c(rep(1, nrow(dat1)), rep(0, nrow(dat0)))
          unlist(calc_mcc(HMeasure(real_labels, c(dat1[["pred"]], dat0[["pred"]]))[["metrics"]]))
        })
        cbind(enc = 1L:ncol(enc_cv), t(enc_cv))
      })
      
      fold_cv <- do.call(rbind, lapply(1L:length(fold_cv), function(single_fold_id)
        cbind(fold = rep(single_fold_id, nrow(fold_cv[[single_fold_id]])), fold_cv[[single_fold_id]])))
    })
    
    rep_cv <- do.call(rbind, lapply(1L:length(rep_cv), function(single_rep_id)
      cbind(repetition = rep(single_rep_id, nrow(rep_cv[[single_rep_id]])), rep_cv[[single_rep_id]])))
  })
  len_cv2 <- do.call(rbind, lapply(1L:length(len_cv), function(single_len_id)
    cbind(len = rep(single_len_id, nrow(len_cv[[single_len_id]])), len_cv[[single_len_id]])))
  
  len_cv2 <- data.frame(len_cv2)
  len_cv2[["len"]] <- as.factor(len_cv2[["len"]])
  levels(len_cv2[["len"]]) <- c("6", "<11", "<16")
  len_cv2[["repetition"]] <- as.factor(len_cv2[["repetition"]])
  len_cv2[["enc"]] <- as.factor(len_cv2[["enc"]])
  len_cv2[["fold"]] <- as.factor(len_cv2[["fold"]])
  
  
  res <- len_cv2 %>% group_by(len, enc, repetition, fold) %>% summarize(mAUC = mean(AUC), mSens = mean(Sens), 
                                                                        mSpec = mean(Spec), mMcc = mean(mcc)) %>%
    group_by(len, enc) %>% summarize(mAUC = mean(mAUC), mSens = mean(mSens), mSpec = mean(mSpec), mMcc = mean(mMcc)) %>%
    filter(mAUC > quantile(mAUC, 0.95)) 
  res3 <- res[order(res[["mAUC"]], decreasing = TRUE), ]
  
  save(list(res, res2, res3), file = "frame_summarize_functions.RData")
  
