source("read_data.R")

source("aa_encodings.R")

library(biogram)
library(cvTools)
library(randomForest)

pos_data <- all_seqs[as.logical(targets), ]
neg_data <- all_seqs[!as.logical(targets), ]


coded_pos <- lapply(aa_groups[1L:2], function(single_group) {
  bitrigrams <- as.matrix(count_multigrams(seq = degenerate(as.matrix(pos_data), single_group), 
                   ns = c(1, rep(2, 4), rep(3, 3)), a()[-1], 
                   ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0))))
  
  bitrigrams <- bitrigrams > 0
  storage.mode(bitrigrams) <- "integer"
  
  bitrigrams
  })

  
coded_neg <- lapply(aa_groups[1L:2], function(single_group) {
  bitrigrams <- as.matrix(count_multigrams(seq = degenerate(as.matrix(neg_data), single_group), 
                                           ns = c(1, rep(2, 4), rep(3, 3)), a()[-1], 
                                           ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                           u = 1L:length(single_group)))
  
  bitrigrams <- bitrigrams > 0
  storage.mode(bitrigrams) <- "integer"
  
  bitrigrams
})
  



fold_res <- pblapply(1L:50, function(dummy) {
  pos_ids <- cvFolds(sum(targets == 1), K = 5)
  neg_ids <- cvFolds(sum(targets != 1), K = 5)
  lapply(1L:5, function(fold) {
    lapply(1L:length(aa_groups), function(group_id) {
      train_pos <- as.matrix(coded_pos[[group_id]])[pos_ids[[5]] != fold, ]
      train_neg <- as.matrix(coded_neg[[group_id]])[neg_ids[[5]] != fold, ]
      
      test_pos <- as.matrix(coded_pos[[group_id]])[pos_ids[[5]] == fold, ]
      test_neg <- as.matrix(coded_neg[[group_id]])[neg_ids[[5]] == fold, ]
      
      test_bis <- test_features(c(rep(1, sum(pos_ids[[5]] != fold)), rep(0, sum(neg_ids[[5]] != fold))),
                                rbind(train_pos, train_neg))
      imp_bigrams <- cut(test_bis, breaks = c(0, 0.05, 1))[[1]]
      
      
      model_cv <- randomForest(x = rbind(train_pos, train_neg)[, imp_bigrams])
      
      preds <- cbind(predict(model_cv, rbind(test_pos, test_neg)), c(rep(1, sum(pos_ids[[5]] == fold)),
                                                                     rep(0, sum(neg_ids[[5]] == fold))))
      preds
    })
  })
})