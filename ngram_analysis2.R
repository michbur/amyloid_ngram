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


res <- len_cv2 %>% group_by(len, enc, repetition, fold) %>% summarize(mAUC = mean(AUC)) %>%
  group_by(len, enc) %>% summarize(mAUC = mean(mAUC)) %>% ungroup 

#best 5 encodings for each length
res %>% group_by(len) %>% arrange(desc(mAUC)) %>% slice(1L:5)

#ids of best encodings
benc_id <- res %>% group_by(len) %>% arrange(desc(mAUC)) %>% slice(1L:5) %>% ungroup %>% 
  select(enc) %>% unlist %>% as.numeric %>% table

#
ggplot(res, aes(y = enc, x = mAUC, colour = len)) +
  geom_point(size = 3) + cool_theme +
  scale_y_discrete("Encoding id") +
  scale_x_continuous("Mean AUC") + scale_fill_discrete("Training length")

#
above080 <- res %>% group_by(enc) %>% summarize(minAUC = min(mAUC)) %>% filter(minAUC > 0.8) %>% 
  select(enc) %>% unlist
ggplot(res %>% filter(enc %in% above080), aes(y = enc, x = mAUC, colour = len)) +
  geom_point(size = 3) + cool_theme +
  scale_y_discrete("Encoding id") +
  scale_x_continuous("Mean AUC") + scale_fill_discrete("Training length")

#
chosen_groups <- groups_summary[groups_summary[["chosen"]] == TRUE, ]
chosen_traits <- chosen_groups[unique(benc_id) - 2, -c(1, 3)]
nice_names <- c("n", "Hydrophobicity", "Solvent accessibility",
                "Polarity", "Interactivity")
colnames(chosen_traits) <- nice_names

#
chosen_traits_tab <- apply(chosen_traits, 2, function(i) data.frame(table(i)))
lapply(1L:length(chosen_traits_tab), function(i) {
  colnames(chosen_traits_tab[[i]]) <- c(nice_names[i], "Frequency")
  chosen_traits_tab[[i]][[1]] <- rownames(aa_nprop)[chosen_traits_tab[[i]][[1]]]
  chosen_traits_tab[[i]]
})

#
apply(chosen_traits[, -1], 2, function(i)
  melt(aa_nprop[unique(i), ]))


save(res, file = "report3.RData")
