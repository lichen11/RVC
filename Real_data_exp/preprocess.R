load("~/Dropbox/Research/Robust_est/Real_data_exp/Celegans/herm_Graph.Rd")
Ag[which(Ag != 0, arr.ind = T)] = 1
pdf("~/Dropbox/Research/Robust_est/Real_data_exp/Figs/AgapAdj.pdf")
pheatmap(Ag,color = colorRampPalette(c("grey", "navy"))(2),
         cluster_row = FALSE, cluster_col = FALSE, display_numbers=F)
dev.off()

pdf("~/Dropbox/Research/Robust_est/Real_data_exp/Figs/AgapAdj2.pdf")
heatmap(Ag, Rowv=NA, Colv=NA, symm=TRUE, cexRow=0.3, cexCol=0.3, add.expr = abline(h=5, v=2))
dev.off()

getTwoCls = function(Ag, vcols, deleteId){
  new.A = rbind(vcols, Ag)
  new.A = new.A[, -which(new.A[1,] == deleteId)]
  new.A = new.A[2:nrow(new.A), ]
  new.A = cbind(vcols, new.A)
  new.A = new.A[-which(new.A[,1] ==deleteId),]
  new.A = new.A[, 2:ncol(new.A)]
  TwoClsA = new.A
  TwoCls = vcols[-which(vcols==deleteId)]
 return(list(TwoClsA, TwoCls))
}
OneTwo = getTwoCls(Ag, vcols, 3)
OneTwoA = OneTwo[[1]]; OneTwoLab = OneTwo[[2]];

OneThree = getTwoCls(Ag, vcols, 2)
OneThreeA = OneThree[[1]]; OneThreeLab = OneThree[[2]];

TwoThree = getTwoCls(Ag, vcols, 1)
TwoThreeA = TwoThree[[1]]; TwoThreeLab = TwoThree[[2]];

writeMat("~/Dropbox/Research/Robust_est/Real_data_exp/Celegans/12ClsMat.mat", OneTwoA = OneTwoA, OneTwoLab = OneTwoLab)
writeMat("~/Dropbox/Research/Robust_est/Real_data_exp/Celegans/13ClsMat.mat", OneThreeA = OneThreeA, OneThreeLab = OneThreeLab)
writeMat("~/Dropbox/Research/Robust_est/Real_data_exp/Celegans/23ClsMat.mat", TwoThreeA = TwoThreeA, TwoThreeLab = TwoThreeLab)
#sum(apply(TwoClsA+t(TwoClsA) , 1 , sum ) == 0)





