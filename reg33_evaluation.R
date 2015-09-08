library(seqinr)
library(randomForest)
library(dplyr)
library(biogram)


#load comparison of rep33 and AmyLoad to remove doubles
compr <- read.csv2("rep33vsAmyLoad.csv", row.names = NULL, stringsAsFactors = FALSE)

#sequences of amyloids existing both in AmyLoad and rep33
doubles <- filter(compr, AmyLoad != "brak") %>% 
  #mutate(nice_name = sapply(strsplit(name, "_", fixed = TRUE), first)) %>%
  select(AmyLoad) %>% unlist(use.names = FALSE) %>% unique

seq_pos <- read.fasta("gcb_abstract_poster/amyloid_pos_full.fasta", seqtype = "AA")

pos_pasted <- sapply(seq_pos, function(i) paste0(i, collapse = ""))

#which sequences should be removed from the positive data set
pos_remove <- which(pos_pasted %in% doubles)

pos_train <- seq_pos[-pos_remove]
neg_train <- read.fasta("gcb_abstract_poster/amyloid_neg_full.fasta", seqtype = "AA")

pos_classifier <- function(pos_dat, neg_dat) {
  lengths_pos <- lengths(pos_dat)
  
}