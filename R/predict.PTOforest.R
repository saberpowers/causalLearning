#'  Make predictions from a fitted PTO forest model
#'
#' @param object a fitted \code{PTOforest} object
#' @param newx matrix of new covariates for which estimated treatment effects
#'  are desired
#' @param ... ignored
#'
#' @return a vector of predictions corresponding to the rows of \code{newx}
#'
#' @export

predict.PTOforest = function(object, newx, ...) {

  colnames(newx) = colnames(object$x)

  if (object$postprocess) {
    stats::predict(object$postfit, data = newx)$predictions

  } else {
    stats::predict(object$PTOfit1, newx = newx) -
      stats::predict(object$PTOfit0, newx = newx)

  }
}

