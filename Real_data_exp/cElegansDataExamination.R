# C.elegans ASE plot
data.dir = "~/Dropbox/Research/Robust_est/Real_data_exp/Celegans/"

load(paste0(data.dir, "herm_Graph.Rd"))
tmp = svd(Ag)
getElbows(tmp$d, 3) #zhu&G uses 13, 87, 176

Xhat = inverse.rdpg(as.matrix(Ag), 13, G = NULL)
Xhat = inverse.rdpg(as.matrix(Ag), 3, G = NULL)
plot3d(Xhat)


mclustModel = Mclust(Xhat, 2) # without specifying, it is 5 componenets
png(paste0(fig.dir, "cElegansASEPlotWithTwoClasses.png"))
plot(Xhat, 
     xlab = "ASE1", ylab = "ASE2")
dev.off()

ASE2d = inverse.rdpg(as.matrix(Ag), 2, G = NULL)
plot(ASE2d)
mclustModelonASE2d = Mclust(ASE2d)

plot(mclustModelonASE2d)


