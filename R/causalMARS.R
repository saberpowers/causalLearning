#'  Fit a causal MARS model
#'
#' @param x matrix of covariates
#' @param tx vector of treatment indicators (0 or 1)
#' @param y vector of response values
#' @param maxterms maximum number of terms to include in the regression basis
#'  (e.g. \code{maxterms = 11} means intercept + 5 pairs added)
#' @param nquant number of quantiles used in splitting
#' @param degree max number of different predictors that can interact in model
#' @param eps shrinkage factor for new term added
#' @param backstep logical: after building out regression basis, should
#'  backward stepwise selection be used to create a sequence of models, with
#'  the criterion evaluated on a validation set to choose among the sequence?
#' @param x.val optional matrix of validation-set covariates
#'  (only used if \code{backstep = TRUE})
#' @param tx.val optional vector of validation-set treatment indicators
#'  (only used if \code{backstep = TRUE})
#' @param y.val optional vector of validation-set response values
#'  (only used if \code{backstep = TRUE})
#' @param propensity logical: should propensity score stratification be used?
#' @param stratum optional vector giving propensity score stratum for each
#'  observation (only used if \code{propensity = TRUE})
#' @param stratum.val optional vector giving propensity score stratum for each
#'  validation-set observation
#'  (only used if \code{propensity = backstep = TRUE})
#' @param minnum minimum number of observations in each arm of each propensity
#'  score stratum needed to estimate regression coefficients for basis
#'  (only used if \code{propensity = TRUE})
#'
#' @return an object of class \code{causalMARS} with attributes:
#'  \itemize{
#'    \item parent: indices of nodes that are parents at each stage
#'    \item childvar: index of predictor chosen at each forward step
#'    \item childquant: quantile of cutoff chosen at each forward step
#'    \item quant: quantiles of the columns of x
#'    \item active: indices of columns with nonzero norm
#'    \item allvars: list of variables appearing in each term
#'    \item khat: the sequence of terms deleted at each step
#'    \item deltahat: relative change in rss
#'    \item rsstesthat: validation-set rss achieved by each model in sequence
#'    \item setesthat: standard error for rsstesthat
#'    \item tim1: time elapsed during forward stepwise phase
#'    \item tim2: total time elapsed
#'    \item x
#'    \item tx
#'    \item y
#'    \item maxterms
#'    \item eps
#'    \item backstep
#'    \item propensity
#'    \item x.val
#'    \item tx.val
#'    \item y.val
#'    \item stratum
#'    \item stratum.val
#'    \item minnum
#'  }
#'
#' @details
#'  parallel arms mars with backward stepwise BOTH randomized case and
#'  propensity stratum. data structures: model terms (nodes) are numbered
#'  1, 2, ... with 1 representing the intercept. forward stepwise:
#'  modmatrix contains basis functions as model is built up -- two columns are
#'  added at each step. Does not include a column of ones for tidiness,
#'  we always add two terms, even when term added in linear (so that reflected
#'  version is just zero).
#'  backward stepwise: khat is the sequence of terms deleted at each step,
#'  based on deltahat = relative change in rss. rsstesthat is rss over test
#'  (validation) set achieved by each reduced model in sequence- used later for
#'  selecting a member of the sequence. active2 contains indices of columns with
#'  nonzero norm 
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
#'# Estimate causal MARS model
#'fit_cm = causalMARS(x, tx, y)
#'pred_cm = predict(fit_cm, newx = x)
#'
#'# Visualize results
#'plot(tx_effect, pred_cm, main = 'Causal MARS',
#'  xlab = 'True treatment effect', ylab = 'Estimated treatment effect')
#'abline(0, 1, lty = 2)
#'
#' @export

causalMARS = function(x, tx, y,
  maxterms = 11, nquant = 5, degree = ncol(x), eps = 1,
  backstep = FALSE, x.val = NULL, tx.val = NULL, y.val = NULL,
  propensity = FALSE, stratum = rep(1, nrow(x)), stratum.val = NULL,
  minnum = 5) {


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


  BIG = 1e+10
  n = nrow(x)
  p = ncol(x)
  x = scale(x, TRUE, FALSE)  #NOTE
  if(!is.null(x.val)) {
    x.val = scale(x.val, TRUE, FALSE)
  }
  
  # compute quantiles for splitting
  discrete = rep(FALSE, p)
  for (j in 1:p) {
    if (length(table(x[, j])) == 2) 
      discrete[j] = TRUE
  }
  
  probs = seq(0, 1, length = nquant)[-nquant]
  quant = apply(x, 2, stats::quantile, probs)
  
  nquantm = rep(nquant, p)
  if (sum(discrete) > 0) {
    for (j in which(discrete)) {
      nquantm[j] = 2
      quant[, j] = NA
      quant[1, j] = 0
    }
  }
  
  if (propensity) {
    stratum.sizes = table(stratum)/nrow(x)
    nstratum = length(stratum.sizes)
    if (sum(as.numeric(names(stratum.sizes)) != (1:length(stratum.sizes))) != 0) 
      stop("Strata should be numbered 1:k")
    stratum.val.sizes = table(stratum.val)/nrow(x.val)
    nstratum.val = length(stratum.val.sizes)
    if (sum(as.numeric(names(stratum.val.sizes)) != (1:length(stratum.val.sizes))) != 
      0) 
      stop("Stratatest should be numbered 1:k")
  }
  
  modmatrix = matrix(0, nrow = n, ncol = maxterms)
  active = rep(FALSE, maxterms)
  active[1] = TRUE
  modmatrix[, 1] = 1
  
  r = y - y
  a0 = c(mean(y[tx == 0]), mean(y[tx == 1]))
  
  r[tx == 0] = y[tx == 0] - mean(y[tx == 0])
  r[tx == 1] = y[tx == 1] - mean(y[tx == 1])
  
  parent = childvar = childquant = NULL
  allvars = vector("list", maxterms)
  maxscorall = rep(NA, maxterms)
  nterms = 1
  
  
  # forward stepwise
  while (nterms < maxterms) {
    maxscor = -1 * BIG
    act = active[1:nterms]
    num = unlist(lapply(allvars, length))[1:nterms]
    act = act & (num < degree)
    
    for (ii in (1:nterms)[act]) {
      jlist = rep(TRUE, p)
      if (ii > 1) 
        jlist[allvars[[ii]]] = FALSE
      
      
      for (j in which(jlist)) {
        for (k in 1:(nquantm[j] - 1)) {
          # bx1 = modmatrix[, ii] * truncpow(x[, j], quant[k, j], dir = 1) bx2 =
          # modmatrix[, ii] * truncpow(x[, j], quant[k, j], dir = 2) NOTE have to make
          # correspondng chanegs in propensity section below!
          bx = NULL
          bx1 = modmatrix[, ii] * truncpow(x[, j], quant[k, j], dir = 1)
          if (sum(bx1^2) > 0) 
          bx = cbind(bx, bx1)
          bx2 = modmatrix[, ii] * truncpow(x[, j], quant[k, j], dir = 2)
          if (sum(bx2^2) > 0) 
          bx = cbind(bx, bx2)
          
          if (!propensity) {
          
            scor = -0.5 * BIG
            if (!is.null(bx)) {
              res = myridge(bx, r, int = TRUE)$res
              res0 = myridge(bx[tx == 0, ], r[tx == 0], int = TRUE)$res
              res1 = myridge(bx[tx == 1, ], r[tx == 1], int = TRUE)$res
              scor = sum(res^2) - sum(res0^2) - sum(res1^2)
            }
          }
          
          if (propensity) {
            if (!is.null(bx)) {
              scor = 0
              
              res0 = res1 = rep(0, n)
              for (s in 1:nstratum) {
                
                res = myridge(bx[stratum == s, ], r[stratum == s])$res
                
                if (sum(tx == 0 & stratum == s) >= minnum) {
                res0[stratum == s & tx == 0] = myridge(bx[tx == 0 & stratum == 
                  s, ], r[tx == 0 & stratum == s])$res
                }
                if (sum(tx == 1 & stratum == s) >= minnum) {
                res1[stratum == s & tx == 1] = myridge(bx[tx == 1 & stratum == 
                  s, ], r[tx == 1 & stratum == s])$res
                }
                scor = scor + (sum(res^2) - sum(res0^2) - sum(res1^2))
                
              }
              res0 = res0[tx == 0]
              res1 = res1[tx == 1]
            }
          }  # end of propensity loop
          
          
          if (scor > maxscor) {
          maxscor = scor
          
          iihat = ii
          jhat = j
          khat = k
          res0hat = res0
          res1hat = res1
          }
        }
      }
    }  #end of for loop

    maxscorall[ii] = maxscor
    new1 = modmatrix[, iihat] * truncpow(x[, jhat], quant[khat, jhat], dir = 1)
    new2 = modmatrix[, iihat] * truncpow(x[, jhat], quant[khat, jhat], dir = 2)
    if (!is.matrix(new1)) 
      new1 = matrix(new1, ncol = 1)
    if (!is.matrix(new2)) 
      new2 = matrix(new2, ncol = 1)
    
    
    modmatrix[, nterms + 1] = new1
    modmatrix[, nterms + 2] = new2
    
    
    active[nterms + 1] = TRUE
    active[nterms + 2] = TRUE
    
    
    allvars[[nterms + 1]] = allvars[[nterms + 2]] = c(allvars[[iihat]], jhat)
    nterms = nterms + 2
    
    
    r[tx == 0] = r[tx == 0] * (1 - eps) + eps * res0hat
    r[tx == 1] = r[tx == 1] * (1 - eps) + eps * res1hat
    
    parent = c(parent, iihat)
    childvar = c(childvar, jhat)
    
    childquant = c(childquant, khat)
  }  #end of while loop
  
  
  active = colSums(modmatrix^2) > 0
  deltahat = khat = NA
  out = list(parent = parent, childvar = childvar, childquant = childquant, quant = quant, 
    active = active, eps = eps, allvars = allvars)
  
  tim1 = proc.time()
  # cat('forward done', fill = TRUE)
  
  
  # backward deletion
  rsstesthat = setesthat = NULL
  
  if (backstep) {

    if(is.null(x.val) | is.null(tx.val) | is.null(y.val)) {
      stop('If backstep = TRUE, then x.val, tx.val, y.val must be specified.')
    }

    BIG = 1e+10
    modmatrix = makebx.newmars(out, x, remove.zerocols = FALSE)[, -1]
    modmatrix.val = makebx.newmars(out, x.val, remove.zerocols = FALSE)[, -1]
    ss = colSums(modmatrix^2) > 0
    
    active2 = rep(TRUE, ncol(modmatrix.val))
    active2[!ss] = FALSE
    khat = deltahat = deltatesthat = rsstesthat = setesthat = rep(NA, sum(active2))
    rtest = rep(NA, nrow(x.val))
    
    go = sum(active2) > 0
    
    ii = 0
    while (go) {
      go = FALSE
      
      delta = deltatest = rsstest = rep(BIG, length(active2))
      # train
      if (!propensity) {
        fit0 = myridge(modmatrix[, active2, drop = FALSE][tx == 0, ], y[tx == 
          0], int = TRUE)
        fit1 = myridge(modmatrix[, active2, drop = FALSE][tx == 1, ], y[tx == 
          1], int = TRUE)
        res0 = fit0$res
        res1 = fit1$res
        
        # test
        yhat0 = cbind(1, modmatrix.val[, active2, drop = FALSE][tx.val == 0, 
          , drop = FALSE]) %*% fit0$coef
        yhat1 = cbind(1, modmatrix.val[, active2, drop = FALSE][tx.val == 1, 
          , drop = FALSE]) %*% fit1$coef
        res0test = (y.val[tx.val == 0] - yhat0)
        res1test = (y.val[tx.val == 1] - yhat1)
        rss0testsq = sum(res0test^2)
        rss1testsq = sum(res1test^2)
      }
      
      if (propensity) 
        {
          rss0sq = rss1sq = rss0testsq = rss1testsq = 0
          # initial fit with parallel terms training
          for (s in 1:nstratum) {
            res0 = res1 = 0
            if (sum(tx == 0 & stratum == s) >= minnum) {
              fit0 = myridge(modmatrix[, active2, drop = FALSE][tx == 0 & stratum == 
              s, ], y[tx == 0 & stratum == s], int = TRUE)
              res0 = fit0$res
            }
            if (sum(tx == 1 & stratum == s) >= minnum) {
              fit1 = myridge(modmatrix[, active2, drop = FALSE][tx == 1 & stratum == 
              s, ], y[tx == 1 & stratum == s], int = TRUE)
              res1 = fit1$res
            }
            rss0sq = rss0sq + sum(res0^2)
            rss1sq = rss1sq + sum(res1^2)
            
            
            # test
            
            res0test = res1test = 0
            if (sum(tx.val == 0 & stratum.val == s) >= minnum) {
              yhat0 = cbind(1, modmatrix.val[, active2, drop = FALSE][tx.val == 
              0 & stratum.val == s, ]) %*% fit0$coef
              res0test = (y.val[tx.val == 0 & stratum.val == s] - yhat0)
            }
            if (sum(tx.val == 1 & stratum.val == s) >= minnum) {
              yhat1 = cbind(1, modmatrix.val[, active2, drop = FALSE][tx.val == 
              1 & stratum.val == s, ]) %*% fit1$coef
              res1test = (y.val[tx.val == 1 & stratum.val == s] - yhat1)
            }
            rss0testsq = rss0testsq + sum(res0test^2)
            rss1testsq = rss1testsq + sum(res1test^2)
          }
        }  #end of propensity
      
      
      for (k in which(active2)) {
        # try collapsing a training set term
        act = active2
        act[k] = FALSE
        
        if (!propensity) {
          
          redfit = myridge(modmatrix[, k, drop = FALSE], y, int = TRUE)
          r = redfit$res
          fit00 = myridge(modmatrix[, act, drop = FALSE][tx == 0, ], r[tx == 
          0], int = TRUE)
          fit11 = myridge(modmatrix[, act, drop = FALSE][tx == 1, ], r[tx == 
          1], int = TRUE)
          
          res00 = fit00$res
          res11 = fit11$res
          delta[k] = (sum(res00^2) + sum(res11^2) - sum(res0^2) - sum(res1^2))/(sum(res00^2) + 
          sum(res11^2))
          # try collapsing the same test set term
          rtest = y.val - cbind(1, modmatrix.val[, k, drop = FALSE]) %*% redfit$coef
          yhat00 = cbind(1, modmatrix.val[, act, drop = FALSE][tx.val == 0, 
          ]) %*% fit00$coef
          yhat11 = cbind(1, modmatrix.val[, act, drop = FALSE][tx.val == 1, 
          ]) %*% fit11$coef
          
          res00test = (rtest[tx.val == 0] - yhat00)
          res11test = (rtest[tx.val == 1] - yhat11)
          rss00testsq = sum(res00test^2)
          rss11testsq = sum(res11test^2)
        }
        
        if (propensity) 
          {
          rss00sq = rss11sq = rss00testsq = rss11testsq = 0
          for (s in 1:nstratum) {
            redfit = myridge(modmatrix[, k, drop = FALSE], y, int = TRUE)
            r = redfit$res
            res00 = res11 = res00test = res11test = 0
            if (sum(tx == 0 & stratum == s) >= minnum) {
            fit00 = myridge(modmatrix[, act, drop = FALSE][tx == 0 & stratum == 
              s, ], r[tx == 0 & stratum == s], int = TRUE)
            res00 = fit00$res
            }
            if (sum(tx == 1 & stratum == s) >= minnum) {
            fit11 = myridge(modmatrix[, act, drop = FALSE][tx == 1 & stratum == 
              s, ], r[tx == 1 & stratum == s], int = TRUE)
            res11 = fit11$res
            }
            
            rss00sq = rss00sq + sum(res00^2)
            rss11sq = rss11sq + sum(res11^2)
            
            # try deleting the same test set term
            
            
            if (sum(tx.val == 0 & stratum.val == s) >= minnum) {
            
            yhat00 = cbind(1, modmatrix.val[, act, drop = FALSE][tx.val == 
              0 & stratum.val == s, ]) %*% fit00$coef
            res00test = (y.val[tx.val == 0 & stratum.val == s] - yhat00)
            }
            if (sum(tx.val == 1 & stratum.val == s) >= minnum) {
            yhat11 = cbind(1, modmatrix.val[, act, drop = FALSE][tx.val == 
              1 & stratum.val == s, ]) %*% fit11$coef
            res11test = (y.val[tx.val == 1 & stratum.val == s] - yhat11)
            }
            
            rss00testsq = rss00testsq + sum(res00test^2)
            rss11testsq = rss11testsq + sum(res11test^2)

          }
          delta[k] = (rss00sq + rss11sq - rss0sq - rss1sq)/(rss00sq + rss11sq)
          }  #end of propensity
        
        # deltatest[k] = (sum(res00test^2) + sum(res11test^2) - sum(res0test^2)
        # -sum(res1test^2))/(sum(res00test^2) + sum(res11test^2))
        
        deltatest[k] = (rss00testsq + rss11testsq - rss0testsq - rss1testsq)/(rss00testsq + 
          rss11testsq)
      }  #end of propensity
      
      ii = ii + 1
      # cat(c('ii=', ii), fill = TRUE)
      khat[ii] = which.min(delta)
      # cat(delta,fill=TRUE)
      deltahat[ii] = delta[khat[ii]]
      deltatesthat[ii] = deltatest[khat[ii]]
      n1 = length(res0test)
      n2 = length(res1test)
      # rsstesthat[ii]=(sum(res0test^2) + sum(res1test^2))/(n1+n2)
      rsstesthat[ii] = (rss0testsq + rss1testsq)/(n1 + n2)
      setesthat[ii] = NA
      # setesthat[ii]=sqrt((n1*var(res0test^2)+n2*var(res1test^2))/(n1+n2)^2) # not
      # currently computed
      
      active2[khat] = FALSE
      num.colsleft = sum(colSums(modmatrix.val[, active2, drop = FALSE]^2) > 0)
      go = (sum(active2) > 0) & (num.colsleft > 0)
      # cat(c(sum(active2), num.colsleft),fill=TRUE) browser() cat(c('dropping ', khat),
      # fill = TRUE)
    }
    
  }
  out$khat = khat[!is.na(khat)]
  out$deltahat = deltahat[!is.na(deltahat)]
  out$rsstesthat = rsstesthat
  out$setesthat = setesthat
  out$tim1 = tim1
  out$tim2 = proc.time()
  out$x = x
  out$tx = tx
  out$y = y
  out$x.val = x.val
  out$tx.val = tx.val
  out$y.val = y.val
  out$maxterms = maxterms
  out$backstep = backstep
  out$propensity = propensity
  out$stratum = stratum
  out$stratum.val = stratum.val
  out$minnum = minnum
  class(out) = 'causalMARS'
  return(out)
}

