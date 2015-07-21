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

len_cv <- pblapply(fold_res, function(single_length) {
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


res <- len_cv2 %>% group_by(len, enc, repetition, fold) %>% 
  summarize(mAUC = mean(AUC), mSpec = mean(Spec), mSens = mean(Sens)) %>%
  group_by(len, enc) %>% summarize(mAUC = mean(mAUC), mSpec = mean(mSpec), mSens = mean(mSens)) %>% ungroup %>% 
  cbind(., n = as.factor(rep(c(4, 4, groups_summary[groups_summary[["chosen"]], "n"]), 3)))

eds <- pbsapply(aa_groups, function(i)
  sapply(aa_groups, function(j)
    calc_ed(j, i)
  )
)

which_len <- "6"

generate_agd <- function(which_len) {
  meds <- melt(eds)
  colnames(meds) <- c("enc1", "enc2", "ed")
  meds[["enc1"]] <- as.factor(meds[["enc1"]])
  meds[["enc2"]] <- as.factor(meds[["enc2"]])
  meds[["mAUC"]] <- filter(res, len == which_len) %>% select(mAUC) %>% unlist %>% rep(95) #%>% cut(breaks = c(0.5, 0.7, 0.8, 0.85))
  meds[["mSpec"]] <- filter(res, len == which_len) %>% select(mSpec) %>% unlist %>% rep(95)
  meds[["mSens"]] <- filter(res, len == which_len) %>% select(mSens) %>% unlist %>% rep(95)
  
  # ggplot(meds, aes(x = enc1, y = enc2, fill = ed, shape = mAUC)) +
  #   geom_tile() +
  #   scale_fill_gradient2(mid = "yellow", high = "red") +
  #   geom_point()
  # looks bad
  
  #aggregated distance data
  agd <- do.call(rbind, lapply(as.character(1L:95), function(i)
    filter(meds, enc2 == i) %>% group_by(ed) %>% summarize(moAUC = mean(mAUC), sdoAUC = sd(mAUC), count = length(mAUC)) %>%
      cbind(enc = rep(i, nrow(.)), .)))
  agd[["ed"]] <- as.factor(agd[["ed"]])
  agd[["count"]] <- cut(agd[["count"]], breaks = c(0, 1, 10, 20, 34))
  #   agd[["moAUC"]] <- cut(agd[["moAUC"]], breaks = c(floor(min(agd[["moAUC"]]) * 100)/100, 
  #                                                    0.7, 0.8, 
  #                                                    ceiling(max(agd[["moAUC"]]) * 100)/100))
  agd[["moAUC"]] <- cut(agd[["moAUC"]], breaks = c(0.6, 
                                                   0.7, 0.8, 
                                                   0.85))
  #agd[["enc"]] <- factor(agd[["enc"]], levels(agd[["enc"]])[order(filter(res, len == which_len)[["mAUC"]], decreasing = TRUE)])
  cbind(len = rep(which_len, nrow(agd)), agd)
}

agds <- do.call(rbind, lapply(c("6", "<11", "<16"), generate_agd))

meds <- melt(eds)
colnames(meds) <- c("enc1", "enc2", "ed")
meds[["enc1"]] <- as.factor(meds[["enc1"]])
meds[["enc2"]] <- as.factor(meds[["enc2"]])




ggplot(agds, aes(x = enc, y = ed, fill = moAUC, shape = count)) +
  geom_tile() +
  geom_point() +
  scale_shape_manual(values = c(3, 4, 8, 15)) +
  facet_wrap(~len, ncol = 1)

ggplot(meds, aes(x = enc1, y = enc2, fill = ed)) +
  geom_tile() +
  scale_fill_gradient2(mid = "yellow", high = "red")

#spec-sens plots
specs <- res %>% group_by(len) %>% filter(mSpec == max(mSpec)) %>% ungroup %>% select(len, enc)
senss <- res %>% group_by(len) %>% filter(mSens == max(mSens)) %>% ungroup %>% select(len, enc)

ss_dat <- cbind(do.call(rbind, lapply(specs[["enc"]], function(i)
  filter(meds, enc1 == i) %>% mutate(enc = enc2, speced = ed) %>% select(enc, speced))),
  do.call(rbind, lapply(senss[["enc"]], function(i)
    filter(meds, enc1 == i) %>% mutate(enc = enc2, sensed = ed) %>% select(sensed))),
  len = unlist(lapply(specs[["len"]], rep, 95)), mAUC = res[, "mAUC"],
  mSpec = res[, "mSpec"], mSens = res[, "mSens"])

ss_dat <- ss_dat %>% group_by(len) %>% 
  mutate(bestAUC = ifelse(mAUC == max(mAUC), "Best AUC", ""),
         bestAUCs = ifelse(mAUC > quantile(mAUC, 0.9), ">0.9 AUC", "<=0.9 AUC")) %>% ungroup

#auc plots
AUCs <- res %>% group_by(len) %>% filter(mAUC == max(mAUC)) %>% ungroup %>% select(len, enc)

AUC_dat <- cbind(do.call(rbind, lapply(AUCs[["enc"]], function(i)
  filter(meds, enc1 == i) %>% mutate(enc = enc2, AUCed = ed) %>% select(enc, AUCed))),
  mAUC = res[, "mAUC"],
  n = as.factor(res[, "n"]), len = res[, "len"])
AUC_dat[["AUCed"]] <- as.factor(AUC_dat[["AUCed"]])

AUC_dat %>% group_by(len) %>% 
  summarise(cor = cor(as.numeric(as.character(AUCed)), mAUC),
            p = cor.test(as.numeric(as.character(AUCed)), mAUC)[["p.value"]])


save(res, meds, agds, ss_dat, AUC_dat, file = "report4.RData")
