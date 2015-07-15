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


