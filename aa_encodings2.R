library(seqinr)
library(dplyr)

data("aaindex")


aa_props <- sapply(aaindex, function(i) i[["I"]])

aa_nprops <- apply(aa_props, 2, function(i) {
  res <- i - min(i, na.rm = TRUE)
  res/max(res, na.rm = TRUE)
})

prop_MK <- read.csv2("AA_index_mk2.csv") %>% filter(!is.na(chosen))

years <- prop_MK %>% select(name) %>% unlist %>% as.character %>% sapply(function(i) 
  strsplit(last(strsplit(i, ", ")[[1]]), ")", fixed = TRUE)[[1]][1])

prop_MK <- cbind(prop_MK, years = years)

prop_MK %>% filter(property == "hydrophobicity") %>% arrange(desc(years))
#c("494", "529", "528")

traits <- prop_MK %>% filter(property != "hydrophobicity") %>% select(X) %>%
  unlist %>% c(., c("494", "529", "528")) %>% as.numeric

best_props <- prop_MK %>% filter(X %in% traits)

grouping_properties <- t(aa_nprop[unlist(traits), ])

all_traits_combn_list <- list(expand.grid(traits))

aa_groups <- unlist(unlist(lapply(all_traits_combn_list, function(all_traits_combn)
  lapply(3L:6, function(single_k)
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
