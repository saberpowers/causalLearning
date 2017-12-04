#'  Make predictions from a bag of fitted causal MARS models
#'
#' @param object a fitted \code{bagged.causalMARS} object
#' @param newx matrix of new covariates for which estimated treatment effects
#'  are desired
#' @param type type of prediction required:
#'  'average' returns a vector of the averages of the bootstrap estimates.
#'  'all' returns a matrix of all of the bootstrap estimates.
#' @param ... ignored
#'
#' @return a vector of estimated personalized treatment effects corresponding
#'  to the rows of \code{newx}
#'
#' @export

predict.bagged.causalMARS = function(object, newx, type = c('average', 'all'),
  ...) {
  type = match.arg(type)
  newx = scale(newx, center = TRUE, scale = FALSE)
  pred = sapply(object, FUN = stats::predict, newx = newx)
  switch(type,
    all = pred,
    average = rowMeans(pred))
}

