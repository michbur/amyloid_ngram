bitrigrams <- as.matrix(count_multigrams(seq = as.matrix(rf_dat[, colnames(rf_dat) != "tar"]), 
                                         ns = c(1, rep(2, 4), rep(3, 3)), a()[-1], ds = list(0, 0, 1, 2, 3, c(0, 0),
                                                                                          c(0, 1), c(1, 0))))
bitrigrams <- bitrigrams > 0
storage.mode(bitrigrams) <- "integer"

test_bis <- test_features(as.numeric(rf_dat[, 11]) - 1, bitrigrams[, -1L:20])
imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]

l1 <- makeLearner("classif.randomForest", id = "bitrigrams", predict.type = "prob")

task_u <- makeClassifTask(id = "bitrigrams", data = data.frame(cbind(bitrigrams[, 1L:20], bitrigrams[, imp_bigrams]), tar = rf_dat[, "tar"]),
                          target = "tar", positive = "pos") 
results <- benchmark(learner = l1, task = task_u, 
                     resampling = makeResampleDesc(method = "CV", iters = 5L, stratify = TRUE), 
                     #auc, sensitivity, specificity
                     measures = list(auc, tnr, tpr),
                     show.info = TRUE)