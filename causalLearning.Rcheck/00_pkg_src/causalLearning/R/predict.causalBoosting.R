#'  Make predictions from a fitted causal boosting model
#'
#' @param object a fitted \code{causalBoosting} object
#' @param newx matrix of new covariates for which estimated treatment effects
#'  are desired
#' @param newtx option vector of new treatment assignments
#'  (only used if \code{type = 'response'})
#' @param type type of prediction required:
#'  'treatment.effect' returns estimated treatment effect.
#'  'conditional.mean' returns two predictions, one for each arm.
#'  'response' returns prediction for arm corresponding to newtx.
#' @param num.trees number(s) of shallow causal trees to use for prediction
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

predict.causalBoosting = function(object, newx, newtx = NULL,
  type = c('treatment.effect', 'conditional.mean', 'response'),
  num.trees = 1:object$num.trees, honest = FALSE, naVal = 0, ...) {
  
  type = match.arg(type)
  if (type == 'response' & is.null(newtx)) {
    stop('response predictions require that newtx be specified')
  }

  CBM = object$CBM

  n = nrow(newx)
  eps = CBM$eps
  G0 = matrix(CBM$intercept[1], n, max(num.trees) + 1)
  G1 = matrix(CBM$intercept[2], n, max(num.trees) + 1)
  for (k in 1:max(num.trees)) {
    pred = stats::predict(CBM$trees[[k]], newx, type = 'conditional.mean',
      honest = honest)
    G1[, k + 1] = G1[, k] + eps * pred[, 2]
    G0[, k + 1] = G0[, k] + eps * pred[, 1]
  }
  G1[is.na(G1)] = naVal
  G0[is.na(G0)] = naVal

  switch(type,
    treatment.effect = (G1 - G0)[, num.trees + 1],
    conditional.mean = list(G1 = G1[, num.trees + 1], G0 = G0[, num.trees + 1]),
    response = (newtx * G1 + (1 - newtx) * G0)[, num.trees + 1])
}

