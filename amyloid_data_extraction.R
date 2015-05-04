#read amyloid data from publication 
dat <- read.csv("amyloid_data.csv")
#remove last row (counts) and middle columns (nothing)
dat <- dat[-nrow(dat), -2]

raw_seqs <- lapply(1L:2, function(single_column) {
  res <- as.character(dat[, single_column])
  res <- res[res != ""]
  #as.vector to unname
  as.vector(sapply(res, function(i)
    substr(i, 2, nchar(i) - 1)))
})


split_raw <- strsplit(c(raw_seqs[[1]], raw_seqs[[2]]), "")


max_len <- max(sapply(split_raw, length))

all_seqs <- t(sapply(split_raw, function(i)
  c(i, rep(NA, 10 - length(i)))))
