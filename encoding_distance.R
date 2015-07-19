source("read_data.R")

source("aa_encodings2.R")

library(biogram)
library(pbapply)

res <- pbsapply(1L:length(aa_groups), function(i)
  sapply(aa_groups[-i], function(j) {
    calc_ed(j, aa_groups[[i]])
  })
)


