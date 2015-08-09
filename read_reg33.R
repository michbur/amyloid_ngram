#read reg33 file

library(dplyr)
library(seqinr)


r33_raw <- readLines("amylpred.33regions.txt") %>% matrix(ncol = 3, byrow = TRUE) %>% 
  data.frame(stringsAsFactors = FALSE)
amyld_ids <- lapply(strsplit(r33_raw[, 3], " "), function(i) matrix(as.numeric(i), nrow = 2))

whole_seqs <- strsplit(r33_raw[[2]], "")

seq_list <- lapply(1L:length(amyld_ids), function(i)
  as.vector(apply(amyld_ids[[i]], 2, function(j) whole_seqs[[i]][j[1]:j[2] + 1])))

all_names <- unlist(lapply(1L:length(r33_raw[[1]]), function(i) {
  names_list <- rep(r33_raw[[1]][i], ncol(amyld_ids[[i]]))
  if(ncol(amyld_ids[[i]]) > 1) {
    paste0(names_list, "_", 1L:ncol(amyld_ids[[i]]))
  } else {
    names_list
  }
}))


seqs <- strsplit(unlist(lapply(seq_list, function(i)
  if(class(i) == "list") {
    sapply(i, paste0, collapse = "")
  } else {
    paste0(i, collapse = "")
  })), "")
write.fasta(seqs, all_names, "reg33.fasta")
