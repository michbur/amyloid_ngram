source("read_data.R")

source("aa_encodings2.R")

library(biogram)

sapply(1L:length(aa_groups), function(i)
  sapply(aa_groups[-i], function(j) {
    print(" ")
    print("LOOP")
    print(" ")
    print(dput(j))
    print(" ")
    print(dput(aa_groups[[i]]))
    calc_ed(j, aa_groups[[i]])
  })
)


a <- structure(list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"), 
                    `2` = c("k", "r", "h"), 
                    `3` = c("d", "e"), 
                    `4` = c("s", "t", "c", "n", "q", "y", "w")), .Names = c("1", "2", "3", "4"))

b <- structure(list(`3` = c("h", "m", "y"),
                    `2` = c("c", "f", "i", "l", "v", "w"), 
                    `1` = c("a", "d", "e", "g", "k", "n", "p", "q", "r", "s", "t")), .Names = c("3", "2", "1"))
calc_ed(a, b)
