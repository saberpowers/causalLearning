#'  Make predictions from a fitted causal MARS model
#'
#' @param object a fitted \code{causalMARS} object
#' @param newx matrix of new covariates for which estimated treatment effects
#'  are desired
#' @param active indices of columns with nonzero norm (defaults to model
#'  selected via backward stepwise phase, or the full model if
#'  \code{backstep = FALSE})
#' @param ... ignored
#'
#' @return a vector of estimated personalized treatment effects corresponding
#'  to the rows of \code{newx}
#'
#' @export

predict.causalMARS = function(object, newx, active, ...) {

  if (missing(active)) {
    if (object$backstep) {
      ntermshat = length(object$khat) - which.min(object$rsstesthat) + 1
      active = (rev(object$khat))[1:ntermshat]
    } else {
      active = 1:(object$maxterms - 1)
    }
  }

  object$x = scale(object$x, center = TRUE, scale = FALSE)
  newx = scale(newx, center = TRUE, scale = FALSE)
  
  del = rep(mean(object$y[object$tx == 1]) - mean(object$y[object$tx == 0]),
    nrow(newx))
  stratum = object$stratum
  stratum.val = object$stratum.val
  minnum = object$minnum
  nstrata = length(unique(stratum))
  stratawt = (table(c(stratum, 1:nstrata)) - 1)
  stratawt = stratawt / sum(stratawt)
  
  if (length(active) > 0) {
    if (!object$propensity) {
      # make training matrices
      bx0 = (makebx.newmars(object, object$x)[, -1])[, active, drop = FALSE]
      beta0 = myridge(bx0[object$tx == 0, ], object$y[object$tx == 0])$coef
      beta1 = myridge(bx0[object$tx == 1, ], object$y[object$tx == 1])$coef
      # make test set matrices
      bx <- (makebx.newmars(object, newx)[, -1])[, active, drop = FALSE]
      del = cbind(1, bx) %*% (beta1 - beta0)
    }
    if (object$propensity) {
      # make training matrices
      bx0 = (makebx.newmars(object, object$x)[, -1])[, active, drop = FALSE]
      beta0 = beta1 = matrix(NA, ncol(bx0) + 1, nstrata)
      for (ss in 1:nstrata) {
        if (sum(object$tx == 0 & stratum == ss) > minnum) {
          beta0[, ss] = myridge(bx0[object$tx == 0 & stratum == ss, ],
            object$y[object$tx == 0 & stratum == ss])$coef
          beta1[, ss] = myridge(bx0[object$tx == 1 & stratum == ss, ],
            object$y[object$tx == 1 & stratum == ss])$coef
        }
      }
      # make test set matrices
      bx <- (makebx.newmars(object, newx)[, -1])[, active, drop = FALSE]
      del = rep(0, nrow(bx))
      totwt = 0
      for (ss in 1:nstrata) {
        if (sum(object$tx == 0 & stratum == ss) > minnum) {
          del = del + stratawt[ss] * cbind(1, bx) %*%
            (beta1[, ss] - beta0[, ss])
          totwt = totwt + stratawt[ss]
        }
      }
      del = del / totwt
    }
  }
  return(c(del))
}

