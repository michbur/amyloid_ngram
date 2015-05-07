bigrams <- cbind(as.matrix(count_ngrams(as.matrix(rf_dat[, -11]), 2, a()[-1], d = 0)),
                 as.matrix(count_ngrams(as.matrix(rf_dat[, -11]), 2, a()[-1], d = 1)),
                 as.matrix(count_ngrams(as.matrix(rf_dat[, -11]), 2, a()[-1], d = 2)),
                 as.matrix(count_ngrams(as.matrix(rf_dat[, -11]), 2, a()[-1], d = 3)),
                 as.matrix(count_ngrams(as.matrix(rf_dat[, -11]), 2, a()[-1], d = 4))) > 0
storage.mode(bigrams) <- "integer"

test_bis <- test_features(as.numeric(rf_dat[, 11]) - 1, 
                          bigrams)
imp_bigrams <- cut(test_bis, breaks = c(0, 0.001, 1))[[1]]

bis <- table_ngrams(as.matrix(rf_dat[, -11]), cut(test_bis, breaks = c(0, 0.05, 1))[[1]], target = as.numeric(rf_dat[, 11]) - 1)

lengths_bis <- 9 - apply(rf_dat[, -11], 1, function(i) sum(is.na(i)))
n_bigrams <- sapply(levels(rf_dat[, 11]), function(i) sum(lengths_bis[rf_dat[, 11] == i]))

bis[["target0"]] <- bis[["target0"]]/n_bigrams["neg"]
bis[["target1"]] <- bis[["target1"]]/n_bigrams["pos"]


mbis <- melt(bis)
#levels(mbis[, "ngram"]) <- decode_ngrams(levels(mbis[, "ngram"]))
levels(mbis[, "variable"]) <- c("nie", "tak")
levels(mbis[, "variable"]) <- c("nie", "tak")

save(mbis, file = "bigram_analysis.RData")
