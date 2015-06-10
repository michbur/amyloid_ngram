library(hmeasure)
library(reshape2)
library(seqinr)
library(ggvis)
library(dplyr)
library(magrittr)
library(ggplot2)

source("aa_encodings.R")
load("amyloid_fold_res.RData")

perf_rep <- function(folds)
  do.call(rbind, lapply(1L:length(folds), function(repetition_id) {
    res <- melt(Reduce("+", lapply(folds[[repetition_id]], function(single_fold)
      sapply(single_fold, function(single_group)
        unlist(HMeasure(single_group[, 2], single_group[, 1])[["metrics"]]))))/5)
    colnames(res) <- c("measure", "encoding", "value")
    res[["encoding"]] <- factor(res[["encoding"]])
    cbind(repetition = factor(rep(repetition_id, nrow(res))), res)
  }))

perf6 <- perf_rep(fold_res6)
perf <- perf_rep(fold_res)




# perf6 %>%
#   filter(measure == "AUC") %>%
#   group_by(encoding) %>%
#   summarize(mAUC = mean(value)) %>%
#   filter(mAUC > quantile(mAUC, 0.9)) %>%
#   droplevels %>%
#   ggvis(y = ~encoding, x = ~mAUC) %>%
#   layer_points()
# 
# perf %>%
#   filter(measure == "AUC") %>%
#   group_by(encoding) %>%
#   summarize(mAUC = mean(value)) %>%
#   filter(mAUC > quantile(mAUC, 0.9)) %>%
#   droplevels %>%
#   ggvis(y = ~encoding, x = ~mAUC) %>%
#   layer_points()

get_best <- function(agg_perf)
  agg_perf %>% 
  filter(measure == "AUC") %>%
  group_by(encoding) %>%
  summarize(mAUC = mean(value)) %>%
  filter(mAUC > quantile(mAUC, 0.9)) %>%
  droplevels %>%
  select(encoding) %>%
  unlist %>%
  as.character %>% 
  as.numeric %>%
  extract(aa_groups, .)

best6 <- perf6 %>% get_best
best <- perf %>% get_best

#get neighbours of amino acid in groupings
#aa - amino acid
get_neigh <- function(group_list) 
  do.call(rbind, lapply(tolower(a()[-1]), function(single_aa) {
    lapply(group_list, function(single_group)
      setdiff(unlist(single_group[sapply(single_group, function(single_subgroup) 
        single_aa %in% single_subgroup)]), single_aa)) %>%
      unlist %>%
      factor(levels = tolower(a()[-1])) %>%
      table(dnn = "n_aa") %>%
      data.frame %>%
      cbind(aa = rep(single_aa, 20), ., perc = .[["Freq"]]/sum(.[["Freq"]]))
  }))

neigh <- rbind(get_neigh(best) %>% cbind(train = rep("all-mers", nrow(.)), .),
               get_neigh(best6) %>% cbind(train = rep("hexamers", nrow(.)), .),
               get_neigh(aa_groups) %>% cbind(train = rep("background", nrow(.)), .))

#unreadable               
# ggplot(neigh, aes(x = n_aa, y = perc, fill = train)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_wrap(~ aa)
# 
# ggplot(neigh, aes(x = aa, y = perc, fill = train)) +
#   geom_bar(stat = "identity") +
#   coord_polar() +
#   facet_wrap(~ n_aa)

cool_theme <- theme(plot.background=element_rect(fill = "transparent",
                                                 colour = "transparent"),
                    panel.grid.major = element_line(colour="lightgrey", linetype = "dashed"),
                    panel.background = element_rect(fill = "transparent",colour = "black"),
                    legend.background = element_rect(fill="NA"))


ggplot(neigh, aes(x = train, y = perc, fill = train)) +
  geom_bar(stat = "identity") +
  facet_grid(n_aa ~ aa) + 
  scale_x_discrete("Amino acid", labels = c("", "", "")) + 
  scale_fill_discrete("Training set") + 
  scale_y_continuous("Frequency") + 
  cool_theme
  
ggplot(neigh, aes(x = train, y = perc, fill = train)) +
  geom_bar(stat = "identity") +
  facet_grid(n_aa ~ aa) + 
  scale_x_discrete(labels = c("", "", "")) + 
  coord_cartesian(ylim = c(0, 0.5)) +
  cool_theme
