#'  Fit a pollinated transformed outcome (PTO) forest model
#'
#' @param x matrix of covariates
#' @param tx vector of treatment indicators (0 or 1)
#' @param y vector of response values
#' @param pscore vector of propensity scores
#' @param num.trees number of trees for transformed outcome forest
#' @param mtry number of variables to possibly split at in each node
#' @param min.node.size minimum node size for transformed outcome forest
#' @param postprocess logical: should optional post-processing random forest be
#'  fit at end?
#' @param verbose logical: should progress be printed to console?
#'
#' @return an object of class \code{PTOforest} with attributes:
#'  \itemize{
#'    \item x: matrix of covariates supplied by function call
#'    \item pscore: vector of propensity score supplied by function call
#'    \item postprocess: logical supplied by function call
#'    \item TOfit: fitted random forest on transformed outcomes
#'    \item PTOfit1: TOfit pollinated with treatment-arm outcomes
#'    \item PTOfit0: TOfit pollinated with control-arm outcomes
#'    \item postfit: post-processing random forest summarizing results
#'  }
#'
#' @examples
#'# Randomized experiment example
#'
#'n = 100 # number of training-set patients to simulate
#'p = 10  # number of features for each training-set patient
#'
#'# Simulate data
#'x = matrix(rnorm(n * p), nrow = n, ncol = p) # simulate covariate matrix
#'tx_effect = x[, 1] + (x[, 2] > 0) # simple heterogeneous treatment effect
#'tx = rbinom(n, size = 1, p = 0.5) # random treatment assignment
#'y = rowMeans(x) + tx * tx_effect + rnorm(n, sd = 0.001) # simulate response
#'
#'# Estimate PTO forest model
#'fit_pto = PTOforest(x, tx, y)
#'pred_pto = predict(fit_pto, newx = x)
#'
#'# Visualize results
#'plot(tx_effect, pred_pto, main = 'PTO forest',
#'  xlab = 'True treatment effect', ylab = 'Estimated treatment effect')
#'abline(0, 1, lty = 2)
#'
#' @export


PTOforest = function(x, tx, y, pscore = rep(.5, nrow(x)),
  num.trees = 500, mtry = ncol(x), min.node.size = max(25, nrow(x) / 40),
  postprocess = TRUE, verbose = FALSE) {


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


  colnames(x) = paste('x', 1:ncol(x), sep = '')
  fit = list(x = x, pscore = pscore, postprocess = postprocess)

  z = tx * y / pscore - (1 - tx) * y / (1 - pscore)

  
  if (verbose) cat('fitting IPW treatment forest\n')
  
  data = data.frame(y = z, x = x)
  colnames(data) = c('y', colnames(x))
  fit$TOfit = ranger::ranger(data = data, dependent.variable.name = 'y',
    num.trees = num.trees, min.node.size = min.node.size, mtry = mtry,
    write.forest = TRUE)

  
# Now pollinate the tree separately with treated and untreated
  if (verbose) {
    cat('pollinating IPW treatment forest separately with treated and',
      'untreated y\n')
  }
  
  fit$PTOfit1 = pollinated.ranger(fit$TOfit, x = x[tx == 1, ], y = y[tx == 1])
  fit$PTOfit0 = pollinated.ranger(fit$TOfit, x = x[tx == 0, ], y = y[tx == 0])
  

  if (postprocess) {
# and one more summarization rf
    if (verbose) cat('fitting TX summary forest\n')

    delta = stats::predict(fit$PTOfit1, x) - stats::predict(fit$PTOfit0, x)
    data = data.frame(y = delta, x = x)
    colnames(x) = paste('x', 1:ncol(x), sep = '')
    colnames(data) = c('y', colnames(x))
    fit$postfit = ranger::ranger(data = data, dependent.variable.name = 'y',
      num.trees = num.trees, mtry = ncol(x), write.forest = TRUE)
  }

  class(fit) = 'PTOforest'
  fit
}

