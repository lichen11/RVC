# some configurations
fig.dir = "~/Dropbox/Research/Robust_est/Real_data_exp/Figs/"

Aall = read.csv(file="~/Dropbox/Research/Robust_est/Real_data_exp/enron/Aall-184x184(noRowColNames).csv",
                header=F, sep=",") 
names(Aall) = NULL
rownames(Aall) = NULL
colnames(Aall) = NULL

# let's use stfp and mclust to get labels
source("~/Dropbox/Research/Local embedding/stfp.R")
source("~/Dropbox/Research/Robust_est/Real_data_exp/getElbows.R")
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
    X <- svd.extract(A, dim, scaling)
  }
  
  # X.mclust <- Mclust(X, G)
  # return(list(X = X, cluster = X.mclust))
  return(X = X)
}

diagAugmentation <- function(A){
  n = nrow(A)
  s <- rowSums(A)
  L <- diag(s)/(n-1) + A  
  return(L)
}
svd.extract = function(A, dim = NULL,scaling = FALSE){
  
  L <- nonpsd.laplacian(A)
  L.svd <- svd(L)
  if(is.null(dim))
    dim <- scree.thresh(L.svd$d)
  
  L.svd.values <- L.svd$d[1:dim]
  L.svd.vectors <- L.svd$v[,1:dim]
  
  if(scaling == TRUE){
    if(dim == 1)
      L.coords <- sqrt(L.svd.values) * L.svd.vectors
    else
      L.coords <- L.svd.vectors %*% diag(sqrt(L.svd.values))
  }
  else{
    L.coords <- L.svd.vectors
  }
  
  return(L.coords)
}
# symmetrize and zero-diagonalize the adjacency matrix
require("sna")
AallSym = symmetrize(as.matrix(Aall), rule = "lower")
#diag(AallSym) = 0
singVal = svd(AallSym)
getElbows(singVal$d, 3) #37  100 138
######### Note, without setting the diagonal to zero, the first elbow is 10.

Xhat = inverse.rdpg(diagAugmentation(AallSym), 10, G = NULL)
mclustModel = Mclust(Xhat, 2)       # without specifying, it is 5 componenets
png(paste0(fig.dir, "enronASEPlotWithTwoClasses.png"))
plot(Xhat, col = mclustModel$classification+1, pch = mclustModel$classification, 
     xlab = "ASE1", ylab = "ASE2")
dev.off()

ASE2d = inverse.rdpg(as.matrix(Aall), 2, G = NULL)
plot(ASE2d)
mclustModelonASE2d = Mclust(ASE2d)
pdf(paste0(fig.dir, "enronMclustOnASE2d4Clust.pdf"))
plot(mclustModelonASE2d)
dev.off()
# now let's do the dark art by augmentating the diagonals

#plot(mclustModel)

require("pheatmap")
png(paste0(fig.dir, "unsortedAdjOverAllTimes.png"))
pheatmap(Aall, cluster_rows=F, cluster_cols=F)
dev.off()

tau = mclustModel$classification
tauDF = data.frame(1:nrow(AallSym), tau)
tauDF = tauDF[order(tau),]
sortedAallSym = AallSym[tauDF[,1], tauDF[,1]]

diag(sortedAall) = 0
png(paste0(fig.dir, "sortedAdjOverAllTimes.png"))
pheatmap(sortedAallSym, cluster_rows=F, cluster_cols=F, show_rownames = F, color = tauDF[,2]+4)
dev.off()

require("R.matlab")

writeMat("adjAndTauOverAllTimes.mat", sortedAall = sortedAallSym,
         sortedTau = tauDF[,2], sortedVertices = tauDF[,1])


install.packages("igraph")
require("igraph")
plot(graph.adjacency(sortedAall, mode = "undirected", diag =T),
     vertex.color= tauDF[,2] + 4, edge.color = "white")

write.table(mclustModel$classification, file = "~/Dropbox/Research/graphMathcing/test/enronLabelsOnAall.csv")

rankInfo = read.csv("~/Dropbox/Research/Robust_est/Real_data_exp/enron/twoClsLabByMclust.csv", header = F, sep = "")

rankInfo = data.frame(rankInfo)
rankInfo[order(rankInfo[,2]),]

write.table(rankInfo[order(rankInfo[,2]),], file = "~/Dropbox/Research/Robust_est/Real_data_exp/enron/twoClsLabByMclustEnron.csv")


data = read.table("~/Dropbox/Research/graphMathcing/test/enronLabels.csv", header = F, sep = ",")
data= data.frame(data)
data[which(data[,4]==1),][,2]














