# some configurations
fig.dir = "~/Dropbox/Research/Robust_est/Real_data_exp/Figs/"



A1 = read.table(file="~/Dropbox/Research/Robust_est/Real_data_exp/enron/A58-184x184.txt",sep=",")              

A2 = read.table(file="~/Dropbox/Research/Robust_est/Real_data_exp/enron/A132-184x184.txt")
A3 = read.table(file="~/Dropbox/Research/Robust_est/Real_data_exp/enron/A136-184x184.txt")
A4 = read.table(file="~/Dropbox/Research/Robust_est/Real_data_exp/enron/A146-184x184.txt")

A = A1 + A2 + A3 + A4

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
getElbows(singVal$d, 3)
######### Note, without setting the diagonal to zero, the first elbow is 10.

# make labels
bicVec = c()
for (d in 1: 30) {
Xhat = inverse.rdpg(as.matrix(A), d, G = NULL)
mclustModel = Mclust(Xhat)
bicVec[d] = mclustModel$bic
}
plot(bicVec)

tmp = svd(Aall)
getElbows(tmp$d, 3)
Aall = as.matrix(Aall)
Xhat = inverse.rdpg(as.matrix(Aall), 10, G = NULL)
mclustModel = Mclust(Xhat, 2) # without specifying, it is 5 componenets
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
tauDF = data.frame(1:nrow(A), tau)
tauDF = tauDF[order(tau),]
sortedAall = Aall[tauDF[,1], tauDF[,1]]
sortedAall = as.matrix(sortedAall)
diag(sortedAall) = 0
png(paste0(fig.dir, "sortedAdjOverAllTimes.png"))
pheatmap(sortedAall, cluster_rows=F, cluster_cols=F, show_rownames = F, color = tauDF[,2]+4)
dev.off()

require("R.matlab")

writeMat("adjAndTauOverAllTimes.mat", sortedAall = sortedAall,
         sortedTau = tauDF[,2], sortedVertices = tauDF[,1])

install.packages("igraph")
require("igraph")
plot(graph.adjacency(sortedAall, mode = "undirected", diag =T),
     vertex.color= tauDF[,2] + 4, edge.color = "white")








