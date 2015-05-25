require("R.matlab")
data.dir = "~/Dropbox/FamousNetworks/"
dataset.dir = paste0(data.dir, "polblogs/")
setwd(dataset.dir)

graph_data = readMat('polblogs.mat')

class = graph_data$Label + 1
#class = readMat("Label.mat")
#class = class$Label
mycol = rainbow(max(class + 2))
A = graph_data$Adj
plotmemb(A, class, inorder=FALSE, main="Political Books Adjacency Matrix",drawborder=TRUE,lwd=1.,lwdb=5,lcol=mycol)
    
