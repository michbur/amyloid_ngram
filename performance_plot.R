#get performance plot data to make ggplot2 plots
get_perf_plot <- function(perf) {
  if (length(perf@alpha.values)!=0) perf@alpha.values <-
    lapply(perf@alpha.values,
           function(x) { isfin <- is.finite(x);
                         x[is.infinite(x)] <-
                           (max(x[isfin]) +
                              mean(abs(x[isfin][-1] -
                                         x[isfin][-length(x[isfin])])));
                         x } )
  
  for (i in 1:length(perf@x.values)) {
    ind.bool <- (is.finite(perf@x.values[[i]]) &
                   is.finite(perf@y.values[[i]]))
    
    if (length(perf@alpha.values)>0)
      perf@alpha.values[[i]] <- perf@alpha.values[[i]][ind.bool]
    
    perf@x.values[[i]] <- perf@x.values[[i]][ind.bool]
    perf@y.values[[i]] <- perf@y.values[[i]][ind.bool]
  }
  
  perf.sampled <- perf
  
  alpha.values <- rev(seq(min(unlist(perf@alpha.values)),
                          max(unlist(perf@alpha.values)),
                          length=max( sapply(perf@alpha.values, length))))
  for (i in 1:length(perf.sampled@y.values)) {
    perf.sampled@x.values[[i]] <-
      approxfun(perf@alpha.values[[i]],perf@x.values[[i]],
                rule=2, ties=mean)(alpha.values)
    perf.sampled@y.values[[i]] <-
      approxfun(perf@alpha.values[[i]], perf@y.values[[i]],
                rule=2, ties=mean)(alpha.values)
  }
  
  ## compute average curve
  perf.avg <- perf.sampled
  perf.avg@x.values <- list( rowMeans( data.frame( perf.avg@x.values)))
  perf.avg@y.values <- list(rowMeans( data.frame( perf.avg@y.values)))
  
  data.frame(x = perf.avg@x.values[[1]], y = perf.avg@y.values[[1]])
}


library(mlr)
library(checkmate)
library(ROCR)


preds <- getBMRPredictions(results, getBMRTaskIds(results), as.df = FALSE)[[1L]]

assertList(preds, c("Prediction", "ResampleResult"), min.len = 1L)
cargs = list(measure = "tpr", x.measure = "fpr")
perf.args = list()
cargs = insert(cargs, perf.args)
plot_preds <- lapply(preds, function(x) {
  cargs$prediction.obj = asROCRPrediction(x)
  do.call(ROCR::performance, cargs)
})


plot(plot_preds[[1]], avg = "threshold")

plot_dat <- do.call(rbind, lapply(plot_preds, get_perf_plot))

classif_nmes <- sapply(strsplit(rownames(plot_dat), ".", fixed = TRUE), function(i) i[1])

#size of the n-gram
ngram_size <- substr(classif_nmes, 0, 1)

#binary ngrams?
ngram_bin <- grepl("bin", classif_nmes)

#distance
dists <- as.numeric(sapply(classif_nmes, function(i) substr(i, nchar(i), nchar(i))))
dists[is.na(dists)] <- "NA"
dists <- paste0("Przerwa: ", dists)
dists <- factor(dists, levels = c("Przerwa: NA", "Przerwa: 0", "Przerwa: 1", "Przerwa: 2", "Przerwa: 3", "Przerwa: 4"))

final_plot_dat <- cbind(plot_dat, dists = dists, Binaryzacja = ngram_bin, n = ngram_size)
rownames(final_plot_dat) <- NULL

library(ggplot2)
ggplot(final_plot_dat, aes(x = x, y = y, col = Binaryzacja, fill = Binaryzacja)) +
  geom_line() +
  scale_x_continuous("TPR") + 
  scale_y_continuous("TFR") +
  geom_point(size=4, shape=21, alpha = 0.5) +
  facet_wrap(~ dists)
