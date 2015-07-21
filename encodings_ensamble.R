source("read_data.R")
source("aa_encodings2.R")
load("amyloid_fold_res2.RData")

library(hmeasure)
library(reshape2)
library(seqinr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(biogram)
library(pbapply)

# len_cv <- pblapply(fold_res, function(single_length) {
#   rep_cv <- lapply(single_length, function(single_rep) {
#     fold_cv <- lapply(single_rep, function(single_fold) {
#       enc_cv <- sapply(single_fold, function(single_enc) {
#         dat <- data.frame(single_enc[[1]])
#         dat1 <- dat %>% filter(X2 == 1) %>% group_by(X3) %>% summarize(pred = max(X1))
#         dat0 <- dat %>% filter(X2 == 0) %>% group_by(X3) %>% summarize(pred = max(X1))
#         real_labels <- c(rep(1, nrow(dat1)), rep(0, nrow(dat0)))
#         unlist(HMeasure(real_labels, c(dat1[["pred"]], dat0[["pred"]]))[["metrics"]])
#       })
#       cbind(enc = 1L:ncol(enc_cv), t(enc_cv))
#     })
#     
#     fold_cv <- do.call(rbind, lapply(1L:length(fold_cv), function(single_fold_id)
#       cbind(fold = rep(single_fold_id, nrow(fold_cv[[single_fold_id]])), fold_cv[[single_fold_id]])))
#   })
#   
#   rep_cv <- do.call(rbind, lapply(1L:length(rep_cv), function(single_rep_id)
#     cbind(repetition = rep(single_rep_id, nrow(rep_cv[[single_rep_id]])), rep_cv[[single_rep_id]])))
# })
# len_cv2 <- do.call(rbind, lapply(1L:length(len_cv), function(single_len_id)
#   cbind(len = rep(single_len_id, nrow(len_cv[[single_len_id]])), len_cv[[single_len_id]])))
# 
# len_cv2 <- data.frame(len_cv2)
# len_cv2[["len"]] <- as.factor(len_cv2[["len"]])
# levels(len_cv2[["len"]]) <- c("6", "<11", "<16")
# len_cv2[["repetition"]] <- as.factor(len_cv2[["repetition"]])
# len_cv2[["enc"]] <- as.factor(len_cv2[["enc"]])
# len_cv2[["fold"]] <- as.factor(len_cv2[["fold"]])
# 
# 
# res <- len_cv2 %>% group_by(len, enc, repetition, fold) %>% 
#   summarize(mAUC = mean(AUC), mSpec = mean(Spec), mSens = mean(Sens)) %>%
#   group_by(len, enc) %>% summarize(mAUC = mean(mAUC), mSpec = mean(mSpec), mSens = mean(mSens)) %>% ungroup %>% 
#   cbind(., n = as.factor(rep(c(4, 4, groups_summary[groups_summary[["chosen"]], "n"]), 3)))
# 
# bestSpec <- which.max(res$mSpec)
# bestSens <- which.max(res$mSens)

#' from code commented out above we get encoding with best specificity and best sensitivity
#' bestSpec = 45
#' bestSens = 277
#' 277%%95=87

# predictions made by bestSpec encoding
rep_cv <- lapply(fold_res[[1]], function(single_rep) {
  fold_cv <- lapply(single_rep, function(single_fold) {
    single_enc <- single_fold[[45]]
    dat <- data.frame(single_enc[[1]])
    dat1 <- dat %>% filter(X2 == 1) %>% group_by(X3) %>% summarize(pred = max(X1))
    dat0 <- dat %>% filter(X2 == 0) %>% group_by(X3) %>% summarize(pred = max(X1))
    real_labels <- c(rep(1, nrow(dat1)), rep(0, nrow(dat0)))
    cbind(rbind(dat1, dat0), real_labels)
  })
  fold_cv <- do.call(rbind, lapply(1L:length(fold_cv), function(single_fold_id)
    cbind(fold = rep(single_fold_id, nrow(fold_cv[[single_fold_id]])), fold_cv[[single_fold_id]])))
})
rep_cv <- do.call(rbind, lapply(1L:length(rep_cv), function(single_rep_id)
  cbind(repetition = rep(single_rep_id, nrow(rep_cv[[single_rep_id]])), rep_cv[[single_rep_id]])))
# predictions made by bestSens encoding
rep_cv2 <- lapply(fold_res[[3]], function(single_rep) {
  fold_cv <- lapply(single_rep, function(single_fold) {
    single_enc <- single_fold[[87]]
    dat <- data.frame(single_enc[[1]])
    dat1 <- dat %>% filter(X2 == 1) %>% group_by(X3) %>% summarize(pred = max(X1))
    dat0 <- dat %>% filter(X2 == 0) %>% group_by(X3) %>% summarize(pred = max(X1))
    real_labels <- c(rep(1, nrow(dat1)), rep(0, nrow(dat0)))
    cbind(rbind(dat1, dat0), real_labels)
  })
  fold_cv <- do.call(rbind, lapply(1L:length(fold_cv), function(single_fold_id)
    cbind(fold = rep(single_fold_id, nrow(fold_cv[[single_fold_id]])), fold_cv[[single_fold_id]])))
})
rep_cv2 <- do.call(rbind, lapply(1L:length(rep_cv2), function(single_rep_id)
  cbind(repetition = rep(single_rep_id, nrow(rep_cv2[[single_rep_id]])), rep_cv2[[single_rep_id]])))

#' we propose different ways of finding consensus betweeen two classifiers
#' This part requires further research in literature as well as numerical experiments

# simple elementwise mean of two vectors
consensus1 <- function(x, y) {
  rowMeans(cbind(x,y))
}

#' elementwise mean of two vectors
#' First classifier is considered an expert regarding specificity
#' Second classifier is considered an expert regarding sensitivity 
consensus2 <- function(x, y, q1=0.1, q2=0.9){
  p <- 0.9
  q <- 0.9
  z <- rowMeans(cbind(x,y))
  z[x<q1] <- p*x[x<q1]+(1-p)*y[x<q1]
  z[y>q2] <- q*y[y>q2]+(1-q)*x[y>q2]
  z[x<q1 & y>q2] <- rowMeans(cbind(x,y))[x<q1 & y>q2]
  z
}

# we combine results and create ensamble classifiers according to consensus functions
combined <- cbind(rep_cv, rep_cv2)
combined[["consensus1"]] <- consensus1(combined[,4], combined[,9])
combined[["consensus2"]] <- consensus2(combined[,4], combined[,9], 0.2, 0.8)
combined <- combined[,c(1,2,3,5,4,9,11,12)]
colnames(combined) <- c("repetition", "fold", "id", "real_labels", "bestSpec", "bestSens", "consensus1", "consensus2")
summary(combined)

# we are most interested in AUC
res2 <- combined %>% group_by(repetition, fold) %>% 
  mutate(auc1=unlist(HMeasure(real_labels, bestSpec)[["metrics"]])[3],
         auc2=unlist(HMeasure(real_labels, bestSens)[["metrics"]])[3],
         auc3=unlist(HMeasure(real_labels, consensus1)[["metrics"]])[3],
         auc4=unlist(HMeasure(real_labels, consensus2)[["metrics"]])[3]) %>%
  summarize(mAUCbSpec = mean(auc1), mAUCbSens = mean(auc2),
            mAUCcons1 = mean(auc3), mAUCcons2 = mean(auc4))

colMeans(res2)

#' one specific repetition and fold
#' it shows that no matter which measure we pick, ensamble works better than
#' each classifier (encoding) seperately
x <- combined %>% filter(repetition==2, fold==5)
unlist(HMeasure(x$real_labels, x$bestSpec)[["metrics"]])
unlist(HMeasure(x$real_labels, x$bestSens)[["metrics"]])
unlist(HMeasure(x$real_labels, x$consensus1)[["metrics"]])
unlist(HMeasure(x$real_labels, x$consensus2)[["metrics"]])

# save(res2, combined, bestSpec, bestSens, file = "encodings_ensamble.RData")