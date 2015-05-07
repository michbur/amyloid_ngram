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


# plot proposed by Piotr

binary_pairs <- rbind(c(1, 2), cbind(3L:7, 3L:7+5))

pair_dat <- cbind(do.call(rbind, lapply(1L:nrow(binary_pairs), function(i) {
  res <- data.frame(competent[, binary_pairs[i, ]], tar = tar, peptide = 1L:nrow(competent))
  colnames(res) <- c("binary", "normal", "tar", "peptide")
  res
})), dists = c(rep("NA", nrow(competent)), as.vector(sapply(0L:4, rep, nrow(competent)))))

pair_dat[["normal"]] <- as.factor(pair_dat[["normal"]])
levels(pair_dat[["normal"]]) <- c("nie", "tak")

pair_dat[["binary"]] <- as.factor(pair_dat[["binary"]])
levels(pair_dat[["binary"]]) <- c("nie", "tak")

pair_dat[["dists"]] <- factor(pair_dat[["dists"]], levels = c("NA", 0:4))
levels(pair_dat[["dists"]]) <- paste0("Przerwa: ", levels(pair_dat[["dists"]]))

levels(pair_dat[["tar"]]) <- c("nie", "tak")

save(pair_dat, file = "pair_dat.RData")

ggplot(pair_dat, aes(x = normal, y =  binary, col = tar, fill = tar)) +
  geom_point(position = "jitter", size=3, shape=21, alpha = 0.5) + 
  scale_colour_hue("Amyloid", l=70, c=150) +
  scale_fill_hue("Amyloid", l=70, c=150) + 
  scale_x_discrete("Zliczenia n-gramów") +
  scale_y_discrete("Zbinaryzowane n-gramy") + 
  facet_wrap(~ dists)
