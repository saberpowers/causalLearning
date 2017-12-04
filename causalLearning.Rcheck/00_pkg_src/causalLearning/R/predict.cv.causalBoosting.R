#'  Make predictions from a fitted cross-validated causal boosting model
#'
#' @param object a fitted \code{cv.causalBoosting} object
#' @param newx matrix of new covariates for which estimated treatment effects
#'  are desired
#' @param newtx option vector of new treatment assignments
#'  (only used if \code{type = 'individual'})
#' @param type type of prediction required:
#'  'treatment.effect' returns estimated treatment effect.
#'  'conditional.mean' returns two predictions, one for each arm.
#'  'response' returns prediction for arm corresponding to newtx.
#' @param num.trees number of shallow causal trees to use for prediction
#' @param naVal value with which to replace \code{NA} predictions
#' @param ... ignored
#'
#' @return a vector or matrix of predictions corresponding to the rows of
#'  \code{newx}
#'
#' @export

predict.cv.causalBoosting = function(object, newx, newtx = NULL,
  type = c('treatment.effect', 'conditional.mean', 'response'),
  num.trees = object$num.trees.min.effect, naVal = 0, ...) {

  predict.causalBoosting(object = object, newx = newx, newtx = newtx,
    type = type, num.trees = num.trees, naVal = naVal)
}

