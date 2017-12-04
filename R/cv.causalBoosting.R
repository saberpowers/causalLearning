#'  Fit a causal boosting model with cross validation
#'
#' @param x matrix of covariates
#' @param tx vector of treatment indicators (0 or 1)
#' @param y vector of response values
#' @param num.trees number of shallow causal trees to build
#' @param maxleaves maximum number of leaves per causal tree
#' @param eps learning rate
#' @param splitSpread how far apart should the candidate splits be for the
#'  causal trees? (e.g. \code{splitSpread = 0.1}) means we consider 10 quantile
#'  cutpoints as candidates for making split
#' @param type.measure loss to use for cross validation:
#'  'response' returns mean-square error for predicting response in each arm.
#'  'effect' returns MSE for treatment effect using honest over-fit estimation.
#' @param nfolds number of cross validation folds
#' @param foldid vector of fold membership
#' @param propensity logical: should propensity score stratification be used?
#' @param stratum optional vector giving propensity score stratum for each
#'  observation (only used if \code{propensity = TRUE})
#' @param isConstVar logical: for the causal tree splitting criterion
#'  (T-statistc), should it be assumed that the noise variance is the same in
#'  treatment and control arms?
#'
#' @return an object of class \code{cv.causalBoosting} which is an object of
#'  class \code{causalBoosting} with these additional attributes:
#'  \itemize{
#'    \item num.trees.min: number of trees with lowest CV error
#'    \item cvm: vector of mean CV error for each number of trees
#'    \item cvsd: vector of standard errors for mean CV errors
#'  }
#'
#' @export


cv.causalBoosting = function(x, tx, y,
  num.trees = 500, maxleaves = 4, eps = 0.01, splitSpread = 0.1,
  type.measure = c('effect', 'response'), nfolds = 5, foldid = NULL,
  propensity = FALSE, stratum = NULL, isConstVar = TRUE) {
  
  type.measure = match.arg(type.measure)

  if (is.null(foldid)) foldid = sample(rep(1:nfolds, length = nrow(x)))
  nfolds = length(unique(foldid))
  
  fit = list()
  pred.response = matrix(0, nrow(x), num.trees)
  pred.effect = matrix(0, nrow(x), num.trees)
  pred.refit = rep(NA, nrow(x))
  
  for (k in 1:nfolds) {

    x.val = x[foldid == k, , drop = FALSE]
    tx.val = tx[foldid == k]
    y.val = y[foldid == k]
    stratum.val = stratum[foldid == k]

    fit[[k]] = causalBoosting(x = x[foldid != k, , drop = FALSE],
      tx = tx[foldid != k], y = y[foldid != k],
      num.trees = num.trees, maxleaves = maxleaves, eps = eps,
      splitSpread = splitSpread, x.est = x.val, tx.est = tx.val, y.est = y.val,
      propensity = propensity, stratum = stratum[foldid != k],
      stratum.est = stratum.val, isConstVar = isConstVar)

    tmp = stats::predict(fit[[k]], x[foldid == k, , drop = FALSE], 1:num.trees,
      type = 'conditional.mean')

    tmpMat = rbind(tmp$G0, tmp$G1)
    n_out = sum(foldid == k)
    pred.response[foldid == k, ] =
      tmpMat[tx[foldid == k] * n_out + (1:n_out), ]

    pred.effect[foldid == k, ] = tmp$G1 - tmp$G0
    tmp.refit = stats::predict(fit[[k]], x[foldid == k, , drop = FALSE],
      1:num.trees, type = 'conditional.mean', honest = TRUE)
    pred.refit[foldid == k] = (tmp.refit$G1 - tmp.refit$G0)[, num.trees]
  }

  cvm.effect = colMeans((pred.effect - pred.refit)^2)
  cvsd.effect = apply((pred.effect - pred.refit)^2, 2, stats::sd) /
    sqrt(nrow(pred.effect))

  cvm.response = apply(pred.response, 2, function(yhat) mean((yhat - y)^2))
  cvsd.response = apply(pred.response, 2,
    function(yhat) stats::sd((yhat - y)^2)) / sqrt(nrow(pred.response))
 
  fit = causalBoosting(x = x, tx = tx, y = y, num.trees = num.trees,
    maxleaves = maxleaves, eps = eps, splitSpread = splitSpread,
    propensity = propensity, stratum = stratum, isConstVar = isConstVar)

  fit$num.trees.min.effect = which.min(cvm.effect)
  fit$num.trees.min.response = which.min(cvm.response)
  fit$num.trees.1se.effect =
    max(which(cvm.effect < min(cvm.effect + cvsd.effect)))

  fit$cvm.effect = cvm.effect
  fit$cvsd.effect = cvsd.effect
  fit$cvm.response = cvm.response
  fit$cvsd.response = cvsd.response

  class(fit) = 'cv.causalBoosting'
  fit
}

