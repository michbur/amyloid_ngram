#random Forests which binarizes their input

#' @export
makeRLearner.classif.randomForestb = function() {
  makeRLearnerClassif(
    cl = "classif.randomForestb",
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
      makeLogicalLearnerParam(id = "keep.inbag", default = FALSE)
    ),
    properties = c("twoclass", "multiclass", "numerics", "factors", "ordered", "prob"),
    name = "Random Forest",
    short.name = "rf"
  )
}


trainLearner.classif.randomForestb = function(.learner, .task, .subset, .weights = NULL, classwt = NULL, cutoff, ...) {
  f = getTaskFormula(.task)
  levs = .task$task.desc$class.levels
  n = length(levs)
  if (missing(cutoff))
    cutoff = rep(1/n, n)
  if (!missing(classwt) && is.numeric(classwt) && length(classwt) == n && is.null(names(classwt)))
    names(classwt) = levs
  if (is.numeric(cutoff) && length(cutoff) == n && is.null(names(cutoff)))
    names(cutoff) = levs
  
  #what an ugly workaound!
  binary_dat <- getTaskData(.task, .subset)
  tmp_dat <- as.matrix(binary_dat[, colnames(binary_dat) != "tar"]) > 0
  storage.mode(tmp_dat) <- "integer"
  
  randomForest::randomForest(f, data = data.frame(tmp_dat, tar = binary_dat[, "tar"]), classwt = classwt, cutoff = cutoff, ...)
}


predictLearner.classif.randomForestb = function(.learner, .model, .newdata, ...) {
  type = ifelse(.learner$predict.type=="response", "response", "prob")
  
  tmp_dat <- as.matrix(.newdata) > 0
  storage.mode(tmp_dat) <- "integer"
  
  predict(.model$learner.model, newdata = data.frame(tmp_dat), type = type, ...)
}


