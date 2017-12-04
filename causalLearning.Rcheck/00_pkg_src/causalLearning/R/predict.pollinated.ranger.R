#'  Make predictions from a pollinated ranger random forest model
#'
#' @param object a fitted \code{pollinated.ranger} object
#' @param newx matrix of new covariates for which predictions are desired
#' @param predict.all logical: should predictions from all trees be returned?
#'  Otherwise the average across trees is returned
#' @param na.treatment how to treat NA predictions from individual trees:
#'  'omit' only uses trees for which the prediction is not NA.
#'  'replace' replaces NA predictions with the overall mean response.
#'  'NA' returns NA if any tree prediction is NA.
#' @param ... additional arguments passed on to \code{predict.ranger}

#'
#' @return a vector of predicted treatment effects corresponding to the rows of
#'  newx

predict.pollinated.ranger = function(object, newx, predict.all = FALSE,
  na.treatment = c('omit', 'replace', 'NA'), ...) {

  na.treatment = match.arg(na.treatment)

  class(object) = 'ranger'
  preds = stats::predict(object, newx, predict.all = TRUE,
    ...)$predictions
  class(object) = 'pollinated.ranger'

  if (na.treatment == 'replace') {
    wh <- is.na(preds)
    preds[wh] = object$mean
  }

  if (!predict.all) {
    preds = apply(preds, 1, mean, na.rm = (na.treatment == 'omit'))
  }

  preds
}

