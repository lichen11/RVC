# make movies
require("vegan")

movie.dir2 = "~/Dropbox/Research/Robust_est/simulationExp/moviesEmbedTo2D/occTo2D/"

B = matrix(c(0.5, 0.42, 0.42, 0.55), 2,2)
rho = c(0.4, 0.6)
n = 200
tmp = getCentersAndDelta(B, rho)
delta = tmp$delta
centers = tmp$centers
rm(tmp)
theoC1 = getCovMat(centers[1,], rho, delta, centers)
theoC2 = getCovMat(centers[2,], rho, delta, centers)

# ##########eigen slow methods in here I didn't pick the top d in magnitude
#spectralDecomp = eigen(simGraph$adjacency)
#spectralDecompAugmented = eigen(diagAugmentation(simGraph$adjacency))
#ase = spectralDecomp$vectors[,1:2] %*% diag(sqrt(spectralDecomp$values[1:2]))
#aseAugmented = spectralDecompAugmented$vectors[,1:2] %*% diag(sqrt(spectralDecompAugmented$values[1:2]))


p_vec = seq(0 , 1, 0.01)
simGraph = sampleSBMfromTrueLatentPosition(n, rho, centers)
#occDif = matrix(NA, nrow = length(p_vec), ncol = 2)
#spectralDecompProcrustesBase = eigen(simGraph$adjacency)
aseProcrustes = inverse.rdpg(diagAugmentation(simGraph$adjacency), 2)
#aseProcrustes = spectralDecompProcrustesBase$vectors[,1:4] %*% diag(sqrt(spectralDecompProcrustesBase$values[1:4]))  

#par(mfrow = c(2,2))
for (i in 1:length(p_vec)){
  p = p_vec[i]
  adjacency = simGraph$adjacency
  set.seed(88)
  contaminationIndex = sample(1:n, round(n*p))
  adjacency[contaminationIndex, contaminationIndex] = matrix(0, nrow = length(contaminationIndex), ncol = length(contaminationIndex))
  # spectralDecomp = eigen(diagAugmentation(adjacency))
  # spectralDecomp = eigen(adjacency)
  # ase = spectralDecomp$vectors[,1:4] %*% diag(sqrt(spectralDecomp$values[1:4]))  
  ase = inverse.rdpg(diagAugmentation(adjacency), 2)
  alignedASE = procrustes(ase, aseProcrustes)$X.new
  #png(paste0(movie.dir2, sprintf("occ%03d.png", i)), width = 500, height = 500)
  #empCov = getSpectralEstAndCov(adjacency, simGraph$tau, 2)
  plot(alignedASE, col = simGraph$tau + 1, pch = simGraph$tau,  
      cex.lab = 1.5, cex.axis = 2, cex.main =2, main = sprintf("p=%g", p))
  #plot(ase, col = simGraph$tau + 1, pch = simGraph$tau,  
   #    cex.lab = 1.5, cex.axis = 2, cex.main =2, main = sprintf("p=%g", p))
  #dev.off()
}


setwd(movie.dir2)
system("convert -delay 15 *.png occDiagAugASEEmbeddedTo2Dmovie.gif")


# make movie of flip linkage

movie.dir4= "~/Dropbox/Research/Robust_est/simulationExp/moviesEmbedTo2D/flipTo2D/"

for (i in 1:length(p_vec)){
  p = p_vec[i]
  adjacency = simGraph$adjacency
  set.seed(88)
  contaminationIndex = sample(1:n, round(n*p))
  tmp = adjacency[contaminationIndex, contaminationIndex] 
  adjacency[contaminationIndex, contaminationIndex] = 1 - tmp
  #spectralDecomp = eigen(diagAugmentation(adjacency))
  #spectralDecomp = eigen(adjacency)
  #ase = spectralDecomp$vectors[,1:2] %*% diag(sqrt(spectralDecomp$values[1:2]))  
  ase = inverse.rdpg(diagAugmentation(adjacency), 2)
  alignedASE = procrustes(ase, aseProcrustes)$X.new
  
  #empCov = getSpectralEstAndCov(adjacency, simGraph$tau, 2)
  png(paste0(movie.dir4, sprintf("flip%03d.png", i)), width = 500, height = 500)
  plot(alignedASE, col = simGraph$tau + 1, pch = simGraph$tau, xlab = "ASE1", ylab = "ASE2", 
       cex.lab = 1.5, cex.axis = 2, cex.main =2, main = sprintf("p=%g", p))    
  dev.off()
}

setwd(movie.dir4)

system("convert -delay 15 *.png flipTo2DTopEvalInMagnitudemovie.gif")

Xhat1 = inverse.rdpg(simGraph$adjacency, 2)
Xhat2 = inverse.rdpg(adjacency, 2)
res1 = Mclust(Xhat1, 2)
res2 = Mclust(Xhat2, 2)
resmat = matrix(NA, nrow = 200, ncol = 2)
resmat[,1] = res1$z[,1]
resmat[,2] = res2$z[,1]
ensemProb = rowMeans(resmat)


#ensemProb = sort(rowMeans(resmat))
adjustedRandIndex(res2$classification, res1$classification)
res3 = data.frame(1:n, ensemProb)
clas1 = res3[order(-res3[,2]),][1:83,]
simGraph$tau[clas1[,1]]

pred = prediction(res2$z[,1], simGraph$tau)
perf = performance(pred,  'tpr', 'fpr')
plot(perf, lwd = 2, box.lty=7, main = "ROC (Classification Tree)")


abline(a=0,b=1)

