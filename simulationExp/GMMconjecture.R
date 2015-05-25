# scratch test code on Avanti etc's conjecture paper
source("~/Dropbox/Research/Local embedding/stfp.R")
require("Rlab")

#########################################################################################
getCentersAndDelta = function(B, rho) {
  centers = eigen(B)$vector %*% diag(sqrt(eigen(B)$values))
  #centers = diag(sqrt(eigen(B)$values)) %*%  eigen(B)$vector 
  x1 = centers[1,]
  x2 = centers[2,]
  delta = outer(x1,x1)*rho[1] + outer(x2, x2)*rho[2]
  return(list(delta = delta, centers = centers))
}

getCovMat = function(x, rho, delta, centers){
  tmp11 = outer(centers[1,], centers[1,])
  tmpx1 = t(x) %*% centers[1,]
  tmp12 = tmpx1 - (tmpx1)^2
  middle1 = as.numeric(tmp12) * tmp11 

  tmp21 = outer(centers[2,], centers[2,])
  tmpx2 = t(x) %*% centers[2,]
  tmp22 = tmpx2 - (tmpx2)^2
  middle2 = tmp21 * as.numeric(tmp22)
  middle = rho[1]*middle1 + rho[2]*middle2
  theoCov = solve(delta)   %*% middle %*% solve(delta)
  return(theoCov)
}

rg.sample = function(P){
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  set.seed(88)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(A)
}

sampleSBMfromTrueLatentPosition = function(n, rho, centers){
  set.seed(88)
  x1 = centers[1,]
  x2 = centers[2,]
  priorVec  = rbern(n, rho[2])
  tau = priorVec + 1
  trueLatentPosition = priorVec %*% t(x2) + (1 - priorVec) %*% t(x1)
  P = trueLatentPosition %*% t(trueLatentPosition)
  A = rg.sample(P)
  return(list(adjacency = A, tau = tau))
}

getSpectralEstAndCov = function(A, tau, d){
  spectralDecomp = eigen(A)
  ase = spectralDecomp$vectors[,1:d] %*% diag(sqrt(spectralDecomp$values[1:d]))
  c1 = var(sqrt(n)*ase[tau ==1,])
  c2 = var(sqrt(n)*ase[tau ==2,])
  return(list(c1 = c1, c2 = c2))
}

getSpectralNormDiff = function(c1, c2, theoC1, theoC2) {
  dif1 = norm(c1-theoC1, "2")
  dif2 = norm(c2-theoC2, "2")
  return(cbind(dif1, dif2))
}
getSpectralNormDifferenceFullFcn = function(B, rho, n, d){
  tmp = getCentersAndDelta(B, rho)
  delta = tmp$delta
  centers = tmp$centers
  rm(tmp)
  theoC1 = getCovMat(centers[1,], rho, delta, centers)
  theoC2 = getCovMat(centers[2,], rho, delta, centers)
  simGraph = sampleSBMfromTrueLatentPosition(n, rho, centers)
  empCov = getSpectralEstAndCov(simGraph$adjacency, simGraph$tau, d)
  dif1 = norm(empCov$c1-theoC1, "2")
  dif2 = norm(empCov$c2-theoC2, "2")
  return(cbind(dif1, dif2))
}

diagAugmentation <- function(A){
  n = nrow(A)
  s <- rowSums(A)
  L <- diag(s)/(n-1) + A  
  return(L)
}

inverse.rdpg = function(A, dim, G = NULL, scaling = FALSE){
  if(is.list(A)){
    for(i in 1:length(A)){
      if(i == 1){
        X <- svd.extract(A[[i]], dim, scaling)
      }
      else{
        X <- cbind(X, svd.extract(A[[i]], dim, scaling))
      }
    }
  }
  else{
    X <- svd.extract(A,dim, scaling)
  }
  
  #X.mclust <- Mclust(X, G)
  #return(list(X = X, cluster = X.mclust))
  return(X = X)
}
#####################################################################################
B = matrix(c(0.7, 0.32, 0.32, 0.75), 2,2)
rho = c(0.4, 0.6)
n = 200
tmp = getCentersAndDelta(B, rho)
delta = tmp$delta
centers = tmp$centers
rm(tmp)
theoC1 = getCovMat(centers[1,], rho, delta, centers)
theoC2 = getCovMat(centers[2,], rho, delta, centers)
simGraph = sampleSBMfromTrueLatentPosition(n, rho, centers)
#spectralDecomp = eigen(simGraph$adjacency)
#spectralDecompAugmented = eigen(diagAugmentation(simGraph$adjacency))
#ase = spectralDecomp$vectors[,1:2] %*% diag(sqrt(spectralDecomp$values[1:2]))
#ase = inverse.rdpg(simGraph$adjacency,2)
#aseAugmented = spectralDecompAugmented$vectors[,1:2] %*% diag(sqrt(spectralDecompAugmented$values[1:2]))


###########Let's use SVD method (for fastness), and Augmented diagonal
ase = inverse.rdpg(diagAugmentation(simGraph$adjacency), 2)


# increase 4000 to 10000
n_vec = seq(5000, 8000, by = 500)
difNormMat = matrix(NA, nrow = length(n_vec), ncol = 2)
for (i in 1: length(n_vec)){
  n = n_vec[i]
  difNormMat[i,] = getSpectralNormDifferenceFullFcn(B, rho, n, 2)
}
png("~/Dropbox/Research/Robust_est/simulationExp/figs/exampleOfIdealCaseEmbedding.png")
plot(ase, col = simGraph$tau + 1, pch = simGraph$tau, xlab = "ASE1", ylab = "ASE2")
dev.off()

plot(n_vec, difNormMat[,1], type = "l")
lo = loess(difNormMat[,1]~n_vec)

png("~/Dropbox/Research/Robust_est/simulationExp/figs/exampleOfIdealCaseCov1Convergence.png")
plot(n_vec[1:12], predict(lo)[1:12], type = 'l', col = 'red', lwd = 4, xlab = "Number of vertices", cex.lab = 1.5, cex.axis = 2, 
     ylab = "Spectral norm of (EmpCov1 - TheoCov1)")
dev.off()

plot(n_vec[1:16], difNormMat[1:16,2], type = "l")
lo2 = loess(difNormMat[,2]~n_vec)
png("~/Dropbox/Research/Robust_est/simulationExp/figs/exampleOfIdealCaseCov2Convergence.png")
plot(n_vec[1:12], predict(lo2)[1:12], type = 'l', col = 'red', lwd = 4, xlab = "Number of vertices", cex.lab = 1.5, cex.axis = 2, 
     ylab = "Spectral norm of (EmpCov2 - TheoCov2)")
dev.off()
save(n_vec, difNormMat, file = "~/Dropbox/Research/Robust_est/simulationExp/result/difSpectralNorm.RData" )
png("~/Dropbox/Research/Robust_est/simulationExp/figs/exampleOfIdealCaseCovConvergenceCombinedPlot.png")
par(mfrow=c(1,2))
plot(n_vec[1:12], predict(lo)[1:12], type = 'l', col = 'red', lwd = 4, xlab = "Number of vertices", 
     ylab = "Spectral norm of (EmpCov1 - TheoCov1)")
plot(n_vec[1:12], predict(lo2)[1:12], type = 'l', col = 'red', lwd = 4, xlab = "Number of vertices", 
     ylab = "Spectral norm of (EmpCov2 - TheoCov2)")
dev.off()

###############       Classification        For the ideal case
# we will try both lda and qda
set.seed(88)
trainSize = round(2/3 * n)
I = sample(1:n, trainSize)
B = matrix(c(0.7, 0.32, 0.32, 0.75), 2,2)
rho = c(0.4, 0.6)

tmp = getCentersAndDelta(B, rho)
delta = tmp$delta
centers = tmp$centers
rm(tmp)
simGraph = sampleSBMfromTrueLatentPosition(n, rho, centers)
spectralDecomp = eigen(simGraph$adjacency)
ase = spectralDecomp$vectors[,1:2] %*% diag(sqrt(spectralDecomp$values[1:2]))
g1 = lda(simGraph$tau ~., data.frame(ase), subset = (1:n)[-I]) 
result = predict(g1, data.frame(ase[I,]))
Y.pred = result$posterior
Y.predClass = result$class
Y.true = simGraph$tau[I]
acc = mean(Y.predClass == Y.true)

g2 = qda(simGraph$tau ~., data.frame(ase), subset = (1:n)[-I]) 
result = predict(g1, data.frame(ase[I,]))
Y.pred = result$posterior
Y.predClass = result$class
Y.true = simGraph$tau[I]
acc = mean(Y.predClass == Y.true)



# what happens under contamination?
# occlusion
# randomly select p vertices and set the entries in the adj to be 0.


B = matrix(c(0.7, 0.32, 0.32, 0.75), 2,2)
rho = c(0.4, 0.6)
n = 500
tmp = getCentersAndDelta(B, rho)
delta = tmp$delta
centers = tmp$centers
rm(tmp)
theoC1 = getCovMat(centers[1,], rho, delta, centers)
theoC2 = getCovMat(centers[2,], rho, delta, centers)
simGraph = sampleSBMfromTrueLatentPosition(n, rho, centers)
spectralDecomp = eigen(simGraph$adjacency)
ase = spectralDecomp$vectors[,1:2] %*% diag(sqrt(spectralDecomp$values[1:2]))


#p_vec = seq(0 , 0.9, 0.05)
p_vec = c(0.15, 0.47, 0.75, 0.95)
simGraph = sampleSBMfromTrueLatentPosition(n, rho, centers)
#occDif = matrix(NA, nrow = length(p_vec), ncol = 2)
#png("~/Dropbox/Research/Robust_est/simulationExp/figs/occASEVaryingOccRate.png")
par(mfrow = c(2,2))
for (i in 1:length(p_vec)){
  p = p_vec[i]
  adjacency = simGraph$adjacency
  set.seed(88)
  contaminationIndex = sample(1:n, round(n*p))
  adjacency[contaminationIndex, contaminationIndex] = matrix(0, nrow = length(contaminationIndex), ncol = length(contaminationIndex))
  #spectralDecomp = eigen(adjacency)
  #ase = spectralDecomp$vectors[,1:2] %*% diag(sqrt(spectralDecomp$values[1:2]))  
  ase = inverse.rdpg(adjacency, 2)
  #alignedASE = procrustes(ase, aseProcrustes)$X.new
 # empCov = getSpectralEstAndCov(adjacency, simGraph$tau, 2)
  plot(ase, col = simGraph$tau + 1, pch = simGraph$tau, xlab = "ASE1", ylab = "ASE2", 
     cex.lab = 1.5, cex.axis = 2, cex.main =2, main = sprintf("p=%g", p))
#  occDif[i,] = getSpectralNormDiff(empCov$c1, empCov$c2, theoC1, theoC2)
  
}


#dev.off()
lo = loess(occDif[,2]~p_vec)
png("~/Dropbox/Research/Robust_est/simulationExp/figs/occCaseCov2Convergence.png")
plot(p_vec, predict(lo), xlab = "p", ylab = "Spectral norm of (EmpCov2 - TheoCov2)", 
     cex.lab = 1.5, cex.axis = 2, cex.main =2, type = "l", col = "red", lwd = 3)
dev.off()
plot(p_vec, occDif[,1], xlab = "p", ylab = "Spectral norm of (EmpCov1 - TheoCov1)", 
          cex.lab = 1.5, cex.axis = 2, cex.main =2, type = "l")


p_vec = c(0.15, 0.47, 0.75, 0.95)
#p_vec = seq(0 , 0.9, 0.05)
simGraph = sampleSBMfromTrueLatentPosition(n, rho, centers)
flipDif = matrix(NA, nrow = length(p_vec), ncol = 2)
#png("~/Dropbox/Research/Robust_est/simulationExp/figs/flipASEVaryingFlipRate.png")
par(mfrow = c(2,2))
for (i in 1:length(p_vec)){
  p = p_vec[i]
  adjacency = simGraph$adjacency
  set.seed(88)
  contaminationIndex = sample(1:n, round(n*p))
  tmp = adjacency[contaminationIndex, contaminationIndex] 
  adjacency[contaminationIndex, contaminationIndex] = 1 - tmp
  #spectralDecomp = eigen(adjacency)
  #ase = spectralDecomp$vectors[,1:2] %*% diag(sqrt(spectralDecomp$values[1:2]))  
  ase = inverse.rdpg(adjacency, 2) 
  plot(ase, col = simGraph$tau + 1, pch = simGraph$tau, xlab = "ASE1", ylab = "ASE2", cex = 2,
       cex.lab = 1.5, cex.axis = 2, cex.main =2, main = sprintf("p=%g", p))
  #empCov = getSpectralEstAndCov(adjacency, simGraph$tau, 2)
  #flipDif[i,] = getSpectralNormDiff(empCov$c1, empCov$c2, theoC1, theoC2)
}
#dev.off()

lo = loess(flipDif[,2]~p_vec)
png("~/Dropbox/Research/Robust_est/simulationExp/figs/flipCaseCov2Convergence.png")
plot(p_vec, predict(lo), xlab = "p", ylab = "Spectral norm of (EmpCov2 - TheoCov2)", 
     cex.lab = 1.5, cex.axis = 2, cex.main =2, type = "l", col = "red", lwd = 3)
dev.off()
plot(p_vec, occDif[,1], xlab = "p", ylab = "Spectral norm of (EmpCov1 - TheoCov1)", 
     cex.lab = 1.5, cex.axis = 2, cex.main =2, type = "l")














