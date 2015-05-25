# scratch no need to save
B = matrix(c(0.4, 0.5, 0.5, 0.5), 2, 2)

#B is the uncontaminated blockmodel
getContaminationB = function(B){
  k = nrow(B)
  upperHalf = cbind(B,B) 
  lowerHalf = cbind(B, matrix(0, nrow = k, ncol = k))
  lowerHalfFlip = cbind(B, matrix(1, nrow = k, ncol = k) - B)
  occB = rbind(upperHalf, lowerHalf)
  flipB = rbind(upperHalf, lowerHalfFlip)
  return(list(occB = occB, flipB = flipB))
}
result = getContaminationB(B)
rankMatrix(result$flipB)



nonpsd.laplacian <- function(A){
    n = nrow(A)
    s <- rowSums(A)
    L <- diag(s)/(n-1) + A
  
    return(L)
}
