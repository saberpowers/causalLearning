#'  Fit a causal boosting model
#'
#' @useDynLib causalLearning
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
#' @param x.est optional matrix of estimation-set covariates used for honest
#'  re-estimation (ignored if \code{tx.est = NULL} or \code{y.est = NULL})
#' @param tx.est optional vector of estimation-set treatment indicators
#'  (ignored if \code{x.est = NULL} or \code{y.est = NULL})
#' @param y.est optional vector of estimation-set response values
#'  (ignored if \code{x.est = NULL} or \code{y.est = NULL})
#' @param propensity logical: should propensity score stratification be used?
#' @param stratum optional vector giving propensity score stratum for each
#'  observation (only used if \code{propensity = TRUE})
#' @param stratum.est optional vector giving propensity score stratum for each
#'  estimation-set observation (ignored if \code{x.est = NULL} or
#'  \code{tx.est = NULL} or \code{y.est = NULL})
#' @param isConstVar logical: for the causal tree splitting criterion
#'  (T-statistc), should it be assumed that the noise variance is the same in
#'  treatment and control arms?
#'
#' @return an object of class \code{causalBoosting} with attributes:
#'  \itemize{
#'    \item CBM: a list storing the intercept, the causal trees and \code{eps}
#'    \item tauhat: matrix of treatment effects for each patient for each step
#'    \item G1: estimated-treatment conditional mean for each patient
#'    \item G0: estimated-control conditional mean for each patient
#'    \item err.y: training error at each step, in predicting response
#'    \item num.trees: number of trees specified by function call
#'  }
#'
#' @export

causalBoosting = function(x, tx, y, num.trees = 500, maxleaves = 4, eps = 0.01,
  splitSpread = 0.1, x.est = NULL, tx.est = NULL, y.est = NULL,
  propensity = FALSE, stratum = NULL, stratum.est = NULL, 
  isConstVar = TRUE) {

 
  # Input sanitization

  x = as.matrix(x)

  if (nrow(x) != length(tx)) {
    stop('nrow(x) does not match length(tx)')

  } else if (nrow(x) != length(y)) {
    stop('nrow(x) does not match length(y)')

  } else if (!is.numeric(x)) {
    stop('x must be numeric matrix')

  } else if (!is.numeric(y)) {
    stop('y must be numeric (use 0/1 for binary response)')

  } else if (!is.numeric(tx) | length(setdiff(tx, 0:1)) > 0) {
    stop('tx must be vector of 0s and 1s')

  }

 
  # s indices are 0-based
  maxNodes = 2 * maxleaves - 1
  
  # if (usePropensity ^ !is.null(s)) { warnings('Non-consistent options: whether to
  # use propensity score will be based on value of s.') }
  
  if (is.null(stratum)) {
    if (propensity) stop('stratum must be specified if propensity = TRUE')
    stratum = -1
  }
  if (is.null(x.est) || is.null(y.est) || is.null(tx.est)) {
    x.est = y.est = tx.est = stratum.est = -1
    n.est = 1
  } else {
    n.est = nrow(x.est)
  if (is.null(stratum.est)) {
      stratum.est = -1
    }
  }
  
  vtxeff = 0
  
  fit = .C("causalBoosting", as.double(x), as.double(y), as.integer(tx),
    as.double(x.est), as.double(y.est), as.integer(tx.est),
    as.integer(num.trees), as.integer(maxleaves), as.double(eps),
    as.integer(propensity), as.integer(stratum), as.integer(stratum.est),
    as.integer(isConstVar), as.integer(nrow(x)), as.integer(ncol(x)),
    as.integer(n.est), as.double(vtxeff), as.double(splitSpread),
    var = integer(num.trees * maxNodes), val = double(num.trees * maxNodes),
    left = integer(num.trees * maxNodes),
    right = integer(num.trees * maxNodes),
    y0bar = double(1), y1bar = double(1), pred0 = double(num.trees * maxNodes),
    pred1 = double(num.trees * maxNodes), cost = double(num.trees * maxNodes),
    pred0e = double(num.trees * maxNodes),
    pred1e = double(num.trees *  maxNodes), G0 = double(nrow(x)),
    G1 = double(nrow(x)), err.y = double(num.trees), err = double(num.trees),
    tauhat = double(num.trees * nrow(x)), PACKAGE = 'causalLearning')
  
  CBM = list()
  CBM$intercept = c(fit$y0bar, fit$y1bar)
  CBM$trees = list()
  CBM$eps = eps
  
  for (k in 1:num.trees) {
    start = (k - 1) * maxNodes + 1
    end = k * maxNodes
    tree = list(var = fit$var[start:end] + 1, val = fit$val[start:end],
      left = fit$left[start:end] + 1, right = fit$right[start:end] + 1,
      pred0 = fit$pred0[start:end], pred1 = fit$pred1[start:end], 
      cost = fit$cost[start:end], pred0e = fit$pred0e[start:end],
      pred1e = fit$pred1e[start:end])
    class(tree) = "causalTree"
    CBM$trees[[k]] = tree
  }
  result = list(CBM = CBM, tauhat = matrix(fit$tauhat, nrow = nrow(x)),
    G1 = fit$G1, G0 = fit$G0, err.y = fit$err.y, num.trees = num.trees)
  
  class(result) = "causalBoosting"
  
  result
}

