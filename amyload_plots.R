library(seqinr)
library(dplyr)
library(ggplot2)

size_mod <- -4
cool_theme <- theme(plot.background=element_rect(fill = "transparent",
                                                 colour = "transparent"),
                    panel.grid.major = element_line(colour="lightgrey", linetype = "dashed"),
                    panel.background = element_rect(fill = "white",colour = "black"),
                    legend.background = element_rect(fill="NA"),
                    legend.position = "bottom",
                    axis.text = element_text(size=12 + size_mod),
                    axis.title.x = element_text(size=16 + size_mod, vjust = -1), 
                    axis.title.y = element_text(size=16 + size_mod, vjust = 1),
                    strip.text = element_text(size=16 + size_mod, face = "bold", lineheight=0.3),
                    strip.background = element_rect(fill="grey", colour = "black"),
                    legend.text = element_text(size=13 + size_mod), 
                    legend.title = element_text(size=17 + size_mod),
                    plot.title = element_text(size=20 + size_mod))


seq_pos <- read.fasta("gcb_abstract_poster/amyloid_pos_full.fasta", seqtype = "AA")
seq_neg <- read.fasta("gcb_abstract_poster/amyloid_neg_full.fasta", seqtype = "AA")

amyl_db <- rbind(data.frame(len = lengths(seq_pos), amyl = "yes"),
                 data.frame(len = lengths(seq_neg), amyl = "no"))

levels(amyl_db[["amyl"]]) <- c("\nAmyloid\n", "\nNon-amyloid\n")

cairo_pdf("density_wholeAL.pdf", pointsize = 1, width=5.56, height=5.56)
ggplot(amyl_db, aes(x = len)) +
  geom_density(fill = "grey", colour = "black") + 
  facet_wrap(~ amyl) +
  #   scale_fill_manual("", values = c("tan2", "skyblue")) +
  #   scale_colour_manual("", values = c("tan2", "skyblue")) +
  scale_x_continuous("Sequence length\n") +
  scale_y_continuous("Density") +
  cool_theme
dev.off()


cairo_pdf("density_cuttedAL.pdf", pointsize = 1, width=5.56, height=5.56)
ggplot(filter(amyl_db, len < 41), aes(x = len)) +
  geom_density(fill = "grey", colour = "black") + 
  facet_wrap(~ amyl) +
  #   scale_fill_manual("", values = c("tan2", "skyblue")) +
  #   scale_colour_manual("", values = c("tan2", "skyblue")) +
  scale_x_continuous("Sequence length\n") +
  scale_y_continuous("Density") +
  cool_theme
dev.off()

amyl_bp <- mutate(amyl_db, clen = cut(len, c(min(len), 5, 6, 10, 15, max(len)), 
                                      include.lowest = TRUE)) %>% 
  group_by(clen, amyl) %>% 
  summarize(n = length(amyl))

cairo_pdf("barplot_AL.pdf", pointsize = 1, width=5.56, height=5.56)
ggplot(amyl_bp, aes(x = clen, y = n, label = n)) +
  geom_bar(fill = "grey", colour = "black", position = "dodge", stat = "identity") +
  geom_text(vjust = -0.25) +
  facet_wrap(~ amyl) +
  scale_y_continuous("Number of sequences") + 
  scale_x_discrete("Range of sequence length\n") +
  cool_theme #+ theme(legend.position = "right")
dev.off()