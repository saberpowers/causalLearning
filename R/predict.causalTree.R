#'  Make predictions from a fitted causal tree model
#'
#' @param object a fitted \code{causalTree} object
#' @param newx matrix of new covariates for which estimated treatment effects
#'  are desired
#' @param newtx option vector of new treatment assignments
#'  (only used if \code{type = 'response'})
#' @param type type of prediction required:
#'  'treatment.effect' returns estimated treatment effect.
#'  'conditional.mean' returns two predictions, one for each arm.
#'  'response' returns prediction for arm corresponding to newtx.
#' @param honest logical: should honest re-estimates of leaf means be used for
#'  prediction? This requires that \code{x.est, tx.est, y.est} were specified
#'  when the causal boosting model was fit
#' @param naVal value with which to replace \code{NA} predictions
#' @param ... ignored
#'
#' @return a vector or matrix of predictions corresponding to the rows of
#'  \code{newx}
#'
#' @export

predict.causalTree = function(object, newx, newtx = NULL,
  type = c('treatment.effect', 'conditional.mean', 'response'), honest = FALSE,
  naVal = 0, ...) {

  leaf = rep(1, nrow(newx))
  terminal = rep(FALSE, nrow(newx))
  while (sum(!is.na(object$var[leaf])) > 0) {
    var = cbind((1:nrow(newx))[!terminal], object$var[leaf[!terminal]])
    leaf[!terminal] = object$left[leaf[!terminal]] + (newx[var] >= object$val[leaf[!terminal]])
    terminal = is.na(object$var[leaf])
  }
  pred = list()
  pred$leaf = leaf
  if (!honest) {
    pred$pred0 = object$pred0[leaf]
    pred$pred1 = object$pred1[leaf]
  } else {
    pred$pred0 = object$pred0e[leaf]
    pred$pred1 = object$pred1e[leaf]
  }
  type = match.arg(type)
  predMatrix = cbind(pred$pred0, pred$pred1)
  predMatrix[is.na(predMatrix)] = naVal
  out = switch(type,
    treatment.effect = predMatrix[, 2] - predMatrix[, 1],
    conditional.mean = predMatrix,
    response = ifelse(!is.null(newtx),
      predMatrix[cbind(1:nrow(newx), newtx + 1)], rep(NA, nrow(newx))))
  out
}

