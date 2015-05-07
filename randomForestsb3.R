#random Forests which binarizes their input and uses QuiPT ---------------------------------------------------------

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
  
  ng_dat <- as.matrix(count_multigrams(seq = as.matrix(dat[, colnames(dat) != "tar"]), ns = c(2, 3), a()[-1], ds = list(0, 0)))
  
  if(binarize) {
    ng_dat <- ng_dat > 0
    storage.mode(ng_dat) <- "integer"
  }
  
  #awful hack to make predict method work
  binarize <<- binarize
  n_gram <<- n_gram
  distance <<- distance
  
  all_features <- test_features(as.numeric(dat[, "tar"]) - 1, ng_dat)
  
  imp_features <<- cut(all_features, breaks = c(0, 5, 1))[[1]]
  
  dat <- data.frame(ng_dat[, imp_features], tar = dat[, "tar"])
  
  randomForest::randomForest(f, data = dat, classwt = classwt, cutoff = cutoff, ...)
}


predictLearner.classif.randomForest_ngram = function(.learner, .model, .newdata, ...) {
  type = ifelse(.learner$predict.type=="response", "response", "prob")
  
  ng_dat <- as.matrix(count_multigrams(seq = as.matrix(.newdata), ns = c(2, 3), a()[-1], ds = list(0, 0)))
  
  if(binarize) {
    ng_dat <- ng_dat > 0
    storage.mode(ng_dat) <- "integer"
  }
  
  predict(.model$learner.model, newdata = data.frame(ng_dat[, imp_features]), type = type, ...)
}
