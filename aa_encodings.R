load("aa_nprop.RData")

traits <- list(size = c(30, 36, 54),
               hydroph = c(1, 26, 33, 57),
               solvent = c(19, 65, 199))

grouping_properties <- t(aa_nprop[unlist(traits), ])

#second one is without polarity
all_traits_combn_list <- list(expand.grid(traits))

aa_groups <- unlist(unlist(lapply(all_traits_combn_list, function(all_traits_combn)
  lapply(3L:5, function(single_k)
    lapply(1L:nrow(all_traits_combn), function(single_trait_combn) {
      cl <- hclust(dist(t(aa_nprop[unlist(all_traits_combn[single_trait_combn, ]), ])))
      gr <- cutree(cl, k = single_k)
      names(gr) <- tolower(names(gr))
      agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
      names(agg_gr) <- 1L:length(agg_gr)
      agg_gr
    }))), recursive = FALSE), recursive = FALSE)


aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("f", "w", "y", "s", "t", "c", "n", "q"))

aa2 = list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("s", "t", "c", "n", "q", "y", "w"))


aa_groups <- c(list(aa2), list(aa1), aa_groups)
