library(seqinr)
data("aaindex")
aa_props <- sapply(aaindex, function(i) i[["D"]])
write.csv2(data.frame(property = unname(aa_props)), file = "aa_props.csv")
