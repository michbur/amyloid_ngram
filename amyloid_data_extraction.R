#read amyloid data from publication 
dat <- read.csv("amyloid_data.csv")
#remove last row (counts) and middle columns (nothing)
dat <- dat[-nrow(dat), -2]

raw_seqs <- lapply(1L:2, function(single_column) {
  res <- as.character(dat[, single_column])
  res <- res[res != ""]
  #as.vector to unname
  as.vector(sapply(res, function(i)
    substr(i, 2, nchar(i) - 1)))
})


split_raw <- strsplit(c(raw_seqs[[1]], raw_seqs[[2]]), "")


max_len <- max(sapply(split_raw, length))

all_seqs <- t(sapply(split_raw, function(i)
  c(i, rep(NA, 10 - length(i)))))

library(signalHsmm)
library(biogram)
library(seqinr)
library(mlr)

tar <- rep(c("pos", "neg"), sapply(raw_seqs, length))
ng2 <- data.frame(as.matrix(count_ngrams(all_seqs, 1, a()[-1])))
rf_dat <- cbind(ng2, tar = as.factor(tar))

task_u <- makeClassifTask(id = "unigramy", data = rf_dat, target = "tar", positive = "pos")
rf <- resample(learner = makeLearner("classif.randomForest", predict.type = "prob"), task = task_u, 
               resampling = makeResampleDesc(method = "CV", iters = 5L, stratify = TRUE), 
               #auc, sensitivity, specificity
               measures = list(auc, tnr, tpr),
               show.info = TRUE)
