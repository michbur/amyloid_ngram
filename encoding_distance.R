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

a <- structure(list(`2` = c("r", "w", "y"), 
                    `1` = c("a", "g", "p", "s", "t"), 
                    `3` = c("d", "e", "h", "k", "n", "q"), 
                    `4` = c("c", "f", "i", "l", "m", "v")), .Names = c("2", "1", "3", "4"))

b <- structure(list(`2` = c("r", "w", "y"), 
                    `3` = c("c", "f", "i", "l", "m", "v"), 
                    `1` = c("a", "d", "e", "g", "h", "k", "n", "p", "q", "s", "t")), 
               .Names = c("2", "3", "1"))
calc_ed(a, b)
