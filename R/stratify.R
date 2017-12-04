#'  Get propensity strata from propensity scores
#'
#' @param pscore vector of propensity scores
#' @param tx vector of treatment indicators
#' @param min.per.arm minimum number of observations for each arm within each
#'  stratum
#'
#' @return a vector of integers with length equal to the length of pscore,
#'  reporting the propensity stratum corresponding to each propensity score
#'
#' @export

stratify = function(pscore, tx, min.per.arm = 30) {

  stratum = ceiling(10 * pscore)
  cutoffs = sort(unique(stratum/10))
  stratum = as.numeric(as.factor(stratum))

  num.treated = stats::aggregate(tx, list(stratum = stratum), sum)$x
  while(min(num.treated) < min.per.arm & length(unique(stratum)) > 1) {
    stratum1 = which.min(num.treated)
    cutoffs = cutoffs[-which.min(num.treated)]
    neighbors = intersect(stratum, stratum1 + c(-1, 1))
    stratum2 = neighbors[which.min(num.treated[neighbors])]
    stratum[stratum == stratum1] = stratum2
    stratum = as.numeric(as.factor(stratum))
    num.treated = stats::aggregate(tx, list(stratum = stratum), sum)$x
  }

  stratum = 1 + max(stratum) - stratum
  cutoffs = rev(cutoffs)

  num.control = stats::aggregate(1 - tx, list(stratum = stratum), sum)$x
  while(min(num.control) < min.per.arm & length(unique(stratum)) > 1) {
    stratum1 = which.min(num.control)
    cutoffs = cutoffs[-which.min(num.control)]
    neighbors = intersect(stratum, stratum1 + c(-1, 1))
    stratum2 = neighbors[which.min(num.control[neighbors])]
    stratum[stratum == stratum1] = stratum2
    stratum = as.numeric(as.factor(stratum))
    num.control = stats::aggregate(tx, list(stratum = stratum), sum)$x
  }

  cutoffs[1] = 1
  cutoffs = rev(cutoffs)

  list(stratum = 1 + max(stratum) - stratum, cutoffs = cutoffs)
}



