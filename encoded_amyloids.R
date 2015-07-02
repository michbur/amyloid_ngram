library(biogram)

pos_data <- all_seqs[as.logical(targets), ]
neg_data <- all_seqs[!as.logical(targets), ]

pos_6 <- apply(pos_data, 1, function(i) length(i) - sum(is.na(i))) == 6
neg_6 <- apply(neg_data, 1, function(i) length(i) - sum(is.na(i))) == 6


coded_pos <- lapply(aa_groups, function(single_group) {
  bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                           ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                           seq = degenerate(as.matrix(pos_data), single_group),
                                           u = as.character(1L:length(single_group))))
  
  bitrigrams <- bitrigrams > 0
  storage.mode(bitrigrams) <- "integer"
  
  bitrigrams
})


coded_neg <- lapply(aa_groups, function(single_group) {
  bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                           ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                           seq = degenerate(as.matrix(neg_data), single_group),
                                           u = as.character(1L:length(single_group))))
  
  bitrigrams <- bitrigrams > 0
  storage.mode(bitrigrams) <- "integer"
  
  bitrigrams
})

save(coded_pos, coded_neg, file = "encoded_amyloids.RData")