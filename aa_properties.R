library(seqinr)
data("aaindex")
aa_props <- sapply(aaindex, function(i) i[["D"]])


property <- rep("", length(aa_props))

property[grep("hydroph", aa_props, ignore.case = TRUE)] <- "hydrophobicity"
property[grep("heli", aa_props, ignore.case = TRUE)] <- "helix"
property[c(grep("sheet", aa_props, ignore.case = TRUE), 6)] <- "sheet"
property[grep("polar", aa_props, ignore.case = TRUE)] <- "polarity"
property[c(33, 34, 127, 319)] <- "solvent surface area"
property[c(9, 59, 63, 73, 109, 112, 150, 399, 515)] <- "size"
write.csv2(data.frame(name = unname(aa_props), property = property), file = "aa_props.csv")
