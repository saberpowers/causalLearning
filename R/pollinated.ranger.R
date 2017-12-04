#'  Pollinate a fitted ranger random forest model
#'
#' @param object a fitted \code{ranger} object
#' @param x matrix of covariates
#' @param y vector of response values
#'
#' @return an object of class \code{pollinated.ranger} which is a \code{ranger}
#'  object that has been pollinated with the data in (x, y)

pollinated.ranger = function(object, x, y) {

  forest = object$forest
  num.trees = forest$num.trees
  which.list = as.list(seq(num.trees))
  split.values = forest$split.values
  split.varIDs = forest$split.varIDs

  for (i in 1:num.trees) {
    which = match(split.varIDs[[i]], 0, FALSE)
    split.values[[i]][which > 0] = seq(sum(which))
    which.list[[i]] = which
  }

  forest$split.values = split.values
  object$forest = forest
  preds = stats::predict(object, x, predict.all = TRUE)$predictions

  ### Get list of means indexed by the unique terminal node values
  pmeans = apply(preds, 2, function(f, y) tapply(y, f, mean), y)

  ### Now populate these terminal nodes with these values
  for (i in 1:num.trees) {
    which = which.list[[i]]
    repvec = rep(NA, sum(which))
    pmean = pmeans[[i]]
    ids = as.integer(names(pmean))
    repvec[ids] = pmean
    split.values[[i]][which > 0] = repvec
  }

  forest$split.values = split.values
  object$forest = forest
  object$mean = mean(y)
  class(object) = c('pollinated.ranger', 'ranger')
  object
}

