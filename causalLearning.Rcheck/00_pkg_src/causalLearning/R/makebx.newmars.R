
makebx.newmars = function(fit, x, remove.zerocols = FALSE) {
  # creates model matrix from 'fit'
  quant = fit$quant
  parent = fit$parent
  childvar = fit$childvar
  childquant = fit$childquant
  active = fit$active
  nterms = length(parent)
  nterms2 = 1 + 2 * nterms
  
  modmatrix = matrix(0, nrow = nrow(x), ncol = nterms2)
  modmatrix[, 1] = 1
  count = 1
  
  for (ii in 1:nterms) {
    new1 = modmatrix[, parent[ii]] * truncpow(x[, childvar[ii]], quant[childquant[ii], 
      childvar[ii]], dir = 1)
    modmatrix[, count + 1] = new1
    new2 = modmatrix[, parent[ii]] * truncpow(x[, childvar[ii]], quant[childquant[ii], 
      childvar[ii]], dir = 2)
    
    modmatrix[, count + 2] = new2
    count = count + 2
  }
  if (remove.zerocols) 
    modmatrix = modmatrix[, active]
  return(modmatrix)
}

