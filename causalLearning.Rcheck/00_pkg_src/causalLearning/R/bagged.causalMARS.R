#'  Fit a bag of causal MARS models
#'
#' @param x matrix of covariates
#' @param tx vector of treatment indicators (0 or 1)
#' @param y vector of response values
#' @param nbag number of models to bag
#' @param maxterms maximum number of terms to include in the regression basis
#'  (e.g. \code{maxterms = 11} means intercept + 5 pairs added)
#' @param nquant number of quantiles used in splitting
#' @param degree max number of different predictors that can interact in model
#' @param eps shrinkage factor for new term added
#' @param backstep logical: should out-of-bag samples be used to prune each
#'  model? otherwise full regression basis is used for each model
#' @param propensity logical: should propensity score stratification be used?
#' @param stratum optional vector giving propensity score stratum for each
#'  observation (only used if \code{propensity = TRUE})
#' @param minnum minimum number of observations in each arm of each propensity
#'  score stratum needed to estimate regression coefficients for basis
#'  (only used if \code{propensity = TRUE})
#' @param verbose logical: should progress be printed to console?
#'
#' @return an object of class \code{bagged.causalMARS}, which is itself a list
#'  of \code{causalMARS} objects
#'
#' @export

bagged.causalMARS = function(x, tx, y, nbag = 20, maxterms = 11, nquant = 5,
  degree = ncol(x), eps = 1, backstep = FALSE,
  propensity = FALSE, stratum = rep(1, nrow(x)), minnum = 5, verbose = FALSE) {

 
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

 
  x = scale(x, center = TRUE, scale = FALSE)

  fit = list()

  for (b in 1:nbag) {

    if (verbose) cat(c('BAG=', b, '/', nbag), fill = TRUE)

    bag = sample(1:nrow(x), size = nrow(x), replace = TRUE)
    oob = rep(TRUE, nrow(x))
    oob[bag] = FALSE
    fit[[b]] = causalMARS(x = x[bag, ], tx = tx[bag], y = y[bag],
      maxterms = maxterms, nquant = nquant, degree = degree, eps = eps,
      backstep = backstep, x.val = x[oob, ], tx.val = tx[oob], y.val = y[oob],
      propensity = propensity, stratum = stratum[bag],
      stratum.val = stratum[oob], minnum = minnum)
  }
  class(fit) = 'bagged.causalMARS'
  fit
}

