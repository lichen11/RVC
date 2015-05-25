lowerBound = function(tmp1, tmp2, p){
  n = length(tmp1)
  k = round(n*p)
  l = n-k
  fill = replicate(k, 0)
  norm1 = norm(tmp1 - tmp2, '2')
  tmp3 = c(tmp1[1:l], fill)
  tmp4 = c(tmp2[1:l], fill)
  norm2 = norm(tmp3 - tmp4, '2')
  result = (norm2 >  min(p*norm1, (1-p)*norm1))
  return(result)
}
cases = c()
nsim = 50000
res = c()
for (i in 1:nsim){
  x1 = runif(500)
  x2 = runif(500)
  p = runif(1, 0.05, 0.95)
  res[i] = lowerBound(x1, x2, p)
  if (res[i] == F){
    cases = c(cases, x1, x2, p)
  }
}
table(res)