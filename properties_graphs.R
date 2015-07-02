source("read_data.R")

source("aa_encodings2.R")

library(dplyr)

plot_values <- prop_MK %>% select(X) %>% unlist %>% c(., 545L:550) %>% slice(data.frame(aa_nprop), .) %>%
  t 

mplot_values <- melt(plot_values)


colnames(mplot_values) <- c("aa", "id", "value")
mplot_values[["id"]] <- as.factor(mplot_values[["id"]])

add_names <- c("Values of Wc in proteins from class Beta, cutoff 6 A, separation 5 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 8 A, separation 5 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 12 A, separation 5 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 6 A, separation 15 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 8 A, separation 15 (Wozniak-Kotulska, 2014)",
               "Values of Wc in proteins from class Beta, cutoff 12 A, separation 15 (Wozniak-Kotulska, 2014)")

plot_names <- prop_MK %>% select(name) %>% unlist %>% as.character %>% c(., add_names)
plot_id <- 1L:length(plot_names)
names(plot_id) <- plot_names

save(plot_values, mplot_values, plot_id, file = "./property_app/properties.RData")



