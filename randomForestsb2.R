#random Forests which binarizes their input ---------------------------------------------------------

#' @export
makeRLearner.classif.randomForestb = function() {
  makeRLearnerClassif(
    cl = "classif.randomForest_ngram",
    package = "randomForest",
    par.set = makeParamSet(
      makeIntegerLearnerParam(id = "ntree", default = 500L, lower = 1L),
      makeIntegerLearnerParam(id = "mtry", lower = 1L),
      makeLogicalLearnerParam(id = "replace", default = TRUE),
      makeNumericVectorLearnerParam(id = "classwt", lower = 0),
      makeNumericVectorLearnerParam(id = "cutoff", lower = 0, upper = 1),
      makeIntegerLearnerParam(id = "sampsize", lower = 1L),
      makeIntegerLearnerParam(id = "nodesize", default = 1L, lower = 1L),
      makeIntegerLearnerParam(id = "maxnodes", lower = 1L),
      makeLogicalLearnerParam(id = "importance", default = FALSE),
      makeLogicalLearnerParam(id = "localImp", default = FALSE),
      makeLogicalLearnerParam(id = "norm.votes", default = TRUE),
      makeLogicalLearnerParam(id = "keep.inbag", default = FALSE),
      #here goes our parameters
      makeLogicalLearnerParam(id = "binarize", default = FALSE),
      makeIntegerLearnerParam(id = "n_gram", default = 1L),
      makeIntegerLearnerParam(id = "distance", default = 0L)
    ),
    properties = c("twoclass", "multiclass", "numerics", "factors", "ordered", "prob", "missings"),
    name = "Random Forest",
    short.name = "rf"
  )
}


trainLearner.classif.randomForest_ngram = function(.learner, .task, .subset, .weights = NULL, classwt = NULL, 
                                                   cutoff, binarize, n_gram, distance, ...) {
  f = getTaskFormula(.task)
  levs = .task$task.desc$class.levels
  n = length(levs)
  if (missing(cutoff))
    cutoff = rep(1/n, n)
  if (!missing(classwt) && is.numeric(classwt) && length(classwt) == n && is.null(names(classwt)))
    names(classwt) = levs
  if (is.numeric(cutoff) && length(cutoff) == n && is.null(names(cutoff)))
    names(cutoff) = levs
  
  #what an ugly workaound! namin convention: tar for target vector
  dat <- getTaskData(.task, .subset)
  
  ng_dat <- as.matrix(count_ngrams(as.matrix(dat[, colnames(dat) != "tar"]), n_gram, a()[-1], distance))
  
  if(binarize) {
    ng_dat <- ng_dat > 0
    storage.mode(ng_dat) <- "integer"
  }
  dat <- data.frame(ng_dat, tar = dat[, "tar"])
  #awful hack to make predict method work
  binarize <<- binarize
  n_gram <<- n_gram
  distance <<- distance
  
  randomForest::randomForest(f, data = dat, classwt = classwt, cutoff = cutoff, ...)
}


predictLearner.classif.randomForest_ngram = function(.learner, .model, .newdata, ...) {
  type = ifelse(.learner$predict.type=="response", "response", "prob")
  
  ng_dat <- as.matrix(count_ngrams(as.matrix(.newdata), n_gram, a()[-1], distance))
  
  if(binarize) {
    ng_dat <- ng_dat > 0
    storage.mode(ng_dat) <- "integer"
  }
  
  predict(.model$learner.model, newdata = data.frame(ng_dat), type = type, ...)
}

