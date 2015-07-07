library(seqinr)
library(dplyr)

data("aaindex")


aa_props <- sapply(aaindex, function(i) i[["I"]])
tableA <- read.table("tableA.csv", sep = ";", dec = ".", head = TRUE)

aa_nprop <- t(apply(cbind(aa_props, tableA[tableA[["X"]] %>% as.character %>% aaa %>% order, 2L:7]), 2, function(i) {
  res <- i - min(i, na.rm = TRUE)
  res/max(res, na.rm = TRUE)
}))

colnames(aa_nprop) <- seqinr::a(colnames(aa_nprop))

prop_MK <- read.csv2("AA_index_mk2.csv") %>% filter(!is.na(chosen))

years <- prop_MK %>% select(name) %>% unlist %>% as.character %>% sapply(function(i) 
  strsplit(last(strsplit(i, ", ")[[1]]), ")", fixed = TRUE)[[1]][1])

prop_MK <- cbind(prop_MK, years = years)

#prop_MK %>% filter(property == "hydrophobicity") %>% arrange(desc(years))
#c("494", "529", "528")



traits <- list(size = 515,
               hydroph = c(494, 529, 528),
               solvent = c(319, 211, 512),
               polarity = c(22, 321),
               interactivity = 545:546)



grouping_properties <- aa_nprop[unlist(traits), ]

all_traits_combn_list <- list(expand.grid(traits))

aa_groups <- unlist(lapply(all_traits_combn_list, function(all_traits_combn)
  lapply(3L:6, function(single_k)
    lapply(1L:nrow(all_traits_combn), function(single_trait_combn) {
      cl <- hclust(dist(t(aa_nprop[unlist(all_traits_combn[single_trait_combn, ]), ])))
      gr <- cutree(cl, k = single_k)
      names(gr) <- tolower(names(gr))
      agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
      names(agg_gr) <- 1L:length(agg_gr)
      agg_gr
    }))), recursive = FALSE)

#perform for each group length separately
aa_id <- lapply(aa_groups, function(j) {
  sort_gr <- lapply(j, function(i) {
    res <- sapply(i, sort)
    res[order(lengths(res))]
  })
  !duplicated(t(sapply(sort_gr, unlist)))
})

aa_groups <- unlist(lapply(1L:length(aa_id), function(i) {
  aa_groups[[i]][aa_id[[i]]]
}), recursive = FALSE)

aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("f", "w", "y", "s", "t", "c", "n", "q"))

aa2 = list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
           `2` = c("k", "r", "h"), 
           `3` = c("d", "e"), 
           `4` = c("s", "t", "c", "n", "q", "y", "w"))


aa_groups <- c(list(aa2), list(aa1), aa_groups)
