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
source("randomForestsb2.R") #all learners

tar <- rep(c("pos", "neg"), sapply(raw_seqs, length))


# case 1 - 1-grams --------------------------

rf_dat <- data.frame(all_seqs, tar = as.factor(tar))

save(rf_dat, file = "amyloid_seqs.RData")




l1 <- makeLearner("classif.randomForestb", id = "1gram", predict.type = "prob", binarize = FALSE, n_gram = 1L, distance = 0)
#RF, binarized data
l2 <- makeLearner("classif.randomForestb", id = "1gram_bin", predict.type = "prob", binarize = TRUE, n_gram = 1L, distance = 0)

#lerners for 2-grams
learners <- unlist(lapply(c(TRUE, FALSE), function(single_binarize)
  lapply(0L:4, function(distance)
    makeLearner("classif.randomForestb", id = paste0("2gram_", ifelse(single_binarize, "bin", ""), distance), 
                predict.type = "prob", binarize = single_binarize, n_gram = 2L, distance = distance))),
  recursive = FALSE)


task_u <- makeClassifTask(id = "unigrams", data = rf_dat, target = "tar", positive = "pos") 
results <- benchmark(learner = c(list(l1), list(l2), learners), task = task_u, 
                         resampling = makeResampleDesc(method = "CV", iters = 5L, stratify = TRUE), 
                         #auc, sensitivity, specificity
                         measures = list(auc, tnr, tpr),
                         show.info = TRUE)
save(results, file = "results.RData")
