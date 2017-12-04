
truncpow = function(x, cutp, dir = 1) {
  if (dir == 1) 
    out = (x - cutp) * (x > cutp)
  if (dir == 2) 
    out = (cutp - x) * (x < cutp)
  return(out)
}

