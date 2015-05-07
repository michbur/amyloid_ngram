#head(results[[1]][[4]][["pred"]][["data"]])
#areas of competence

load("results.RData")
library(dplyr)

competent <- do.call(cbind, lapply(results[[1]], function(single_learner) 
  single_learner[["pred"]][["data"]] %>% group_by(id) %>% summarise(competent = truth == response) %>% 
    ungroup %>% select(competent)
))

colnames(competent) <- sapply(results[[1]], function(single_learner) single_learner[["learner.id"]])

competent <- as.matrix(competent)
storage.mode(competent) <- "integer"
library(reshape2)

mcompetent <- melt(competent)
ngram_bin <- factor(grepl("bin", as.character(mcompetent[["Var2"]])))


ngram_size <- substr(as.character(mcompetent[["Var2"]]), 0, 1)

dists <- as.numeric(sapply(as.character(mcompetent[["Var2"]]), function(i) substr(i, nchar(i), nchar(i))))
dists[is.na(dists)] <- "NA"
dists <- factor(dists)

#tar is defined in amyloid_data_exraction.R
comp_dat <- cbind(mcompetent, bin = ngram_bin, n = ngram_size, dists = dists, tar = rep(tar, ncol(competent)))

diffs_c <- cbind(abs(filter(comp_dat, bin == FALSE) %>% select(value) - filter(comp_dat, bin == TRUE) %>% select(value)),
      filter(comp_dat, bin == FALSE) %>% select(-value))

aoc <- diffs_c %>% group_by(dists, tar) %>% summarise(count = sum(value)) %>% data.frame
aoc[aoc[["tar"]] == "neg", "count"] <- filter(aoc, tar == "neg") %>% select(count) / sum(tar == "neg")
aoc[aoc[["tar"]] == "pos", "count"] <- filter(aoc, tar == "pos") %>% select(count) / sum(tar == "pos")
aoc[["dists"]] <- factor(aoc[["dists"]], levels = c("NA", 0L:4))
levels(aoc[["tar"]]) <- c("Nie", "Tak")
save(aoc, file = "aoc.RData")


library(ggplot2)

ggplot(aoc, aes(x = dists, y = count, fill = tar)) +
  geom_bar(stat = "identity", position = "dodge")
