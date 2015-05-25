# check scree plots
setwd("~/Dropbox/Research/Robust_est/simulationExp/")
source('~/Dropbox/Research/Robust_est/simulationExp/getElbows.R')
#remember to source some functions in GMMconjection.R
B = matrix(c(0.5, 0.42, 0.42, 0.55), 2,2)
rho = c(0.4, 0.6)
n = 200
tmp = getCentersAndDelta(B, rho)
delta = tmp$delta
centers = tmp$centers
rm(tmp)
theoC1 = getCovMat(centers[1,], rho, delta, centers)
theoC2 = getCovMat(centers[2,], rho, delta, centers)

p_vec = seq(0 , 1, 0.01)
simGraph = sampleSBMfromTrueLatentPosition(n, rho, centers)
aseProcrustes = inverse.rdpg(diagAugmentation(simGraph$adjacency), 2)
scree_selected_mat = matrix(NA, nrow = length(p_vec), ncol = 3)
for (i in 1:length(p_vec)){
    p = p_vec[i]
    adjacency = simGraph$adjacency
    set.seed(i)
    contaminationIndex = sample(1:n, round(n*p))
    adjacency[contaminationIndex, contaminationIndex] = matrix(0, nrow = length(contaminationIndex), ncol = length(contaminationIndex))
    sing_val = svd(adjacency)$d[1:20]
    scree_selected_mat[i,] = getElbows(sing_val, n = 3, threshold = F, plot = F) 
    #ase = inverse.rdpg(diagAugmentation(adjacency), 2)
    #png(paste0(movie.dir2, sprintf("occ%03d.png", i)), width = 500, height = 500)
    #empCov = getSpectralEstAndCov(adjacency, simGraph$tau, 2)
    
    #plot(ase, col = simGraph$tau + 1, pch = simGraph$tau,  
    #    cex.lab = 1.5, cex.axis = 2, cex.main =2, main = sprintf("p=%g", p))
    #dev.off()
}
# fix p
p = 0.74
nsim = 500
simGraph = sampleSBMfromTrueLatentPosition(n, rho, centers)
aseProcrustes = inverse.rdpg(diagAugmentation(simGraph$adjacency), 2)
scree_selected_mat2 = matrix(NA, nrow = nsim, ncol = 3)
for (i in 1:nsim){
    
    adjacency = simGraph$adjacency
    set.seed(i)
    contaminationIndex = sample(1:n, round(n*p))
    adjacency[contaminationIndex, contaminationIndex] = matrix(0, nrow = length(contaminationIndex), ncol = length(contaminationIndex))
    sing_val = svd(adjacency)$d[1:22]
    scree_selected_mat2[i,] = getElbows(sing_val, n = 3, threshold = F, plot = F) 
    #ase = inverse.rdpg(diagAugmentation(adjacency), 2)
    #png(paste0(movie.dir2, sprintf("occ%03d.png", i)), width = 500, height = 500)
    #empCov = getSpectralEstAndCov(adjacency, simGraph$tau, 2)
    
    #plot(ase, col = simGraph$tau + 1, pch = simGraph$tau,  
    #    cex.lab = 1.5, cex.axis = 2, cex.main =2, main = sprintf("p=%g", p))
    #dev.off()
}

