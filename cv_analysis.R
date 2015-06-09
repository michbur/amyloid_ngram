library(hmeasure)

Reduce("+", lapply(fold_res6[[1]], function(single_fold)
  sapply(single_fold, function(single_group)
    unlist(HMeasure(single_group[, 2], single_group[, 1])[["metrics"]]))))/5
