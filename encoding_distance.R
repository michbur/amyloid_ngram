source("read_data.R")

source("aa_encodings2.R")

library(biogram)

sapply(1L:length(aa_groups), function(i)
  sapply(aa_groups[-i], function(j) {
    print(c(j, aa_groups[[i]]))
    calc_ed(j, aa_groups[[i]])
    })
  )

a <- structure(list(`3` = c("c", "i", "l", "v"), 
               `2` = c("f", "h", "m", "r", "w", "y"), 
               `1` = c("a", "d", "e", "g", "k", "n", "p", "q", "s", "t")), .Names = c("3", "2", "1"))

b <- structure(list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
               `2` = c("k", "r", "h"), `3` = c("d", "e"), 
               `4` = c("s", "t", "c", "n", "q", "y", "w")), .Names = c("1", "2", "3", "4"))
calc_ed(b, a)
