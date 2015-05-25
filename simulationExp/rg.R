library("irlba")
library("mclust")
library("cluster")
library("MASS")
library("reshape")
library("plyr")
library("ggplot2")

## Sample an undirected graph on n vertices
## Input: P is n times n matrix giving the parameters of the Bernoulli r.v.
rg.sample <- function(P){
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- runif(n*(n-1)/2)
  U <- (U + t(U))
  A <- (U < P) + 0 ;
  diag(A) <- 0
  return(A)
}

rg.sample.pois <- function(P){
  n <-  nrow(P)
  U <- matrix(0, nrow = n, ncol = n)
  U[col(U) > row(U)] <- rpois(n*(n-1)/2, lambda = P[col(P) > row(P)])
  U <-  U + t(U)
  return(U)
}

rg.SBM <- function(n, B, rho,condition=FALSE){
  if(!condition){
    tau <- sample(c(1:length(rho)), n, replace = TRUE, prob = rho)
  }
  else{
    tau <- unlist(lapply(1:2,function(k) rep(k, rho[k]*n)))
  }
  P <- B[tau,tau]
  return(list(adjacency=rg.sample(P),tau=tau))
}

stfp <- function(A, dim, method = "svd"){

    if(method == "svd"){
        S <- irlba(A, nu = dim, nv = dim)
        V <- S$v[,1:dim]
        D <- S$d[1:dim]
    } else if(method=="eigensmall"){
        n = nrow(A)
        S <- eigen(A)
        V <- S$vectors[,(n-dim+1):n]
        D <- S$values[(n-dim+1):n]      
    }
    else{
        S <- eigen(A)
        V <- S$vectors[,1:dim]
        D <- S$values[1:dim]
    }

    if(dim == 1){
        Xhat <- V*sqrt(D)
    }
    else{
      Xhat <- V %*% diag(sqrt(D))
    }
   
    return(Xhat)
}

procrustes <- function(A,B){
    mm <- svd(t(A) %*% B)
    W <- mm$u %*% t(mm$v)
    return(W)
}

sd.multivariate <- function(X){
    n <- nrow(X)
    Xbar <- outer(rep(1,n), colMeans(X))
    Xtilde <- X - Xbar
    return(1/(n-1)*(t(Xtilde)%*%Xtilde))
}

block.var <- function(X,rho){
    n <-  nrow(X)
    var.list <- list()

    Delta <-  matrix(0, nrow = ncol(X), ncol = ncol(X))
    for( i in 1:n){
        Delta <- Delta + outer(X[i,],X[i,])*rho[i]
    }
   
    for(i in 1:n){
        tmp1 <- X[i,]%*%t(X)
        tmp2 <- tmp1 - tmp1^2
        B <- matrix(0, nrow = ncol(X), ncol = ncol(X))
        for(j in 1:n){
            B <-  B + outer(X[j,], X[j,])*tmp2[j]*rho[j]
        }
        var.list[[i]] <- solve(Delta)%*%B%*%solve(Delta)
    }
    return(var.list)
}

bayes.mvn <- function(n, mu1, mu2, Sigma1, Sigma2, pi1, pi2, nmc){

    tau <- sample(1:2, nmc, replace = TRUE, prob = c(pi1,pi2))
    m1 <- sum(tau == 1)
    m2 <-  sum(tau == 2)

    X1 <- mvrnorm(m1, sqrt(n)*mu1, Sigma1)
    X2 <-  mvrnorm(m2, sqrt(n)*mu2, Sigma2)
    labels <- c(rep(1,m1),rep(-1,m2))

    x <- rbind(X1,X2)

    Sigma1.inv <- solve(Sigma1)
    Sigma2.inv <- solve(Sigma2)

    c1 <- log(pi1) - 1/2*log(det(Sigma1))
    c2 <- log(pi2) - 1/2*log(det(Sigma2))

    xbar1 <- x - sqrt(n)*outer(rep(1,nrow(x)),mu1)
    xbar2 <- x - sqrt(n)*outer(rep(1,nrow(x)),mu2)

    density1 <- c1 - rowSums((1/2*xbar1%*%Sigma1.inv)*xbar1)
    density2 <- c2 - rowSums((1/2*xbar2%*%Sigma2.inv)*xbar2)

    labels.hat <- sign(density1 - density2)
    error <- sum(labels.hat != labels)/nmc
    return(error)
}

   
experiment1 <- function(n){
    #n <- 8000
    rho <- c(0.4,0.6)
    B <- matrix(c(0.7,0.32,0.32,0.75), nrow = 2, ncol = 2)

    A <- rg.SBM(n, B,rho)
    tau <- A$tau
    x <- eigen(B)
    xx <- x$vectors %*% diag(sqrt(x$values))
    X <- xx[tau,]

    Xhat <- stfp(A$adjacency,dim = 2)
    W <- procrustes(Xhat,X)
    residue <- sqrt(n)*(Xhat%*%W - X)
    print(norm(Xhat%*%W - X,"F")^2)

    sigs <- block.var(xx,rho)

    print(sigs)

   #  mean.residue1 <- colMeans(residue[which(tau == 1),])
   sd.residue1 <- sd.multivariate(residue[which(tau == 1),])
   #  mean.residue2 <- colMeans(residue[which(tau == 2),])
   sd.residue2 <- sd.multivariate(residue[which(tau == 2),])
   print(sd.residue1)
   print(sd.residue2)

   # ## pdf(paste("clusplot",n,".pdf",sep=''),useDingbats=FALSE)
   #  aa <- Mclust(sqrt(n)*Xhat,2,modelNames = c("VVV"))
   #  error.rate.gmm <- min(sum(abs(aa$classification - tau)),
   #                    sum(abs(3 - aa$classification - tau)))/n
    df <- data.frame(x = (Xhat%*%W)[,1], y = (Xhat%*%W)[,2],
                     block = as.character(tau),
                     tau = tau)
##    qplot(data = dat, x = X1, y = X2, colour = labels) + stat_ellipse()

    library(ellipse)
    df_ell <- data.frame()

    for(g in unique(df$tau)){
        sig <- sigs[[g]]
        df_ell <- rbind(df_ell, cbind(as.data.frame(
          with(df[df$tau==g,],
               ellipse(sig[1,2]/sqrt(sig[1,1]*sig[2,2]), 
               scale=c(sqrt(sig[1,1])/sqrt(n),sqrt(sig[2,2])/sqrt(n)), 
               centre=c(xx[g,1],xx[g,2]))),group=g)))

    }
    
    library(ggplot2)
    ggplot(data=df, aes(x=x, y=y,color=block)) + 
    geom_point(size=2, alpha=.3) +  theme(panel.background = element_rect(fill = 'white', colour = 'black')) +  
    geom_path(data=df_ell[1:(nrow(df_ell)/2),], aes(x=x, y=y), size=1,linetype = 2, colour=1) +
          geom_path(data = df_ell[-c(1:(nrow(df_ell)/2)),], aes(x = x, y = y), size = 1, linetype = 2, colour = 1)
    ggsave(paste("clusplot",n,".pdf",sep=''),width=6,height=4)

    
    # clusplot(sqrt(n)*Xhat%*%W,aa$classification, color=TRUE, shade = TRUE,
    #          span = FALSE, sub = "",
    #          main = paste("Gaussian mixture clustering of the residue for n = ",
    #            n, sep = ''), xlab = "", ylab = "")
    # dev.off()

    # bb <- kmeans(sqrt(n)*Xhat,2, iter.max = 50)
    # clusplot(sqrt(n)*Xhat%*%W,bb$cluster,color = TRUE, shade = TRUE, sub = "", main = paste("Kmeans clustering of the residue for n = ", n, sep = ''))
    
    # error.rate.kmeans <- min(sum(abs(bb$cluster - tau)),
    #                   sum(abs(3 - bb$cluster - tau)))/n

    return(df)
 
}

experiment2 <- function(){
    set.seed(12345)
    nseq <- seq(from = 1000, to = 4000, by = 250)
    rho <- c(0.6,0.4)
    B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)
    x <- eigen(B)
    xx <- x$vectors %*% diag(sqrt(x$values))

    error.vec.gmm <- numeric(length(nseq))
    error.vec.kmeans <-  numeric(length(nseq))
    error.vec.bayes <- numeric(length(nseq))

    tmp <- block.var(xx, rho)
    Sigma1 <- tmp[[1]]
    Sigma2 <- tmp[[2]]
    nmc <- 100

    for(i in 1:length(nseq)){
        error.i.gmm <- numeric(nmc)
        error.i.kmeans <- numeric(nmc)
        for(j in 1:nmc){
            ## A <- rg.SBM(nseq[i], B, rho)
            tau <- rep(c(1:nrow(B)),nseq[i]*rho)
            P <- B[tau,tau]
            A <- rg.sample(P)
            X <- xx[tau,]

            diag(A) <- rowSums(A)/sqrt(sum(rowSums(A)))

            Xhat <- stfp(A,2,method="eigen")
            aa <- Mclust(sqrt(nseq[i])*Xhat,2,modelNames = c("VVV"))
            tmp <- sum(aa$classification != tau)/nseq[i]
            error.i.gmm[j] <- min(tmp, 1 - tmp)

            bb <- kmeans(sqrt(nseq[i])*Xhat,2, iter.max = 50)
            tmp <- sum(bb$cluster != tau)/nseq[i]
            error.i.kmeans[j] <- min(tmp, 1 - tmp)
        }
        error.vec.gmm[i] <- mean(error.i.gmm)
        error.vec.kmeans[i] <- mean(error.i.kmeans)
        error.vec.bayes[i] <- bayes.mvn(nseq[i], xx[1,], xx[2,], Sigma1, Sigma2, rho[1], rho[2], 100000)
   }
        error.vec.log <- log(nseq)/nseq
        error.vec.log <- error.vec.kmeans[1]/error.vec.log[1]*error.vec.log
 
        dat.new <- melt(data.frame(nseq, gmm = log10(error.vec.gmm), kmeans = log10(error.vec.kmeans), bayes = log10(error.vec.bayes), log = log10(error.vec.log)), id = "nseq")
        ggplot(dat.new, aes(x = nseq, y = value, colour = variable)) + geom_line() + xlab("n") +
          ylab("classification error") + labs(title = element_blank()) +
            theme(panel.background = element_rect(fill = 'white', colour = 'black')) +
              theme(legend.title = element_blank()) +
                theme(legend.key = element_rect(fill = 'white', colour = 'white')) + scale_color_discrete(breaks=c("gmm", "kmeans", "bayes", "log"), labels = c("GMM", "K-Means", "Bayes", "STFP"))
       ggsave("gmm_kmeans_bayes.pdf", width=6,height=4)
    
}
 
# if sum(rho) == 1, randomly allocated block memberships
# if sum(rho) == n, conditionally allocated block memberships
sbm = function(n,K,B,rho,myseed=F){
if(myseed) set.seed(myseed)
if(sum(rho) == n) blockn = rho
if(sum(rho) != n) blockn = tabulate(sample(K,n,replace=T,prob=rho),K)
tau=NULL; for(k in 1:K) tau = c(tau,rep(k,blockn[k]))
A = matrix(0,nrow=n,ncol=n)
for(i in 1:(n-1)) for(j in (i+1):n) A[i,j] = A[j,i] = rbinom(1,1,B[tau[i],tau[j]])
return(list(A,tau))
}

experiment2.dan <- function(nmc=100, nseq=seq(from = 1000, to = 4000, by = 250)){
  rho <- c(0.6,0.4)
  B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)
  x <- eigen(B)
  xx <- x$vectors %*% diag(sqrt(x$values))

  tmp <- block.var(xx, rho)
  Sigma1 <- tmp[[1]]
  Sigma2 <- tmp[[2]]

  
  all <- data.frame(n=numeric(0),error=numeric(0),method=character(0))
  error.max = 0
  for(n in nseq){
    
    for(j in 1:nmc){
      A <- rg.SBM(n, B, rho)
      tau <- A$tau
      X <- xx[tau,]
      
      Xhat <- stfp(A$adjacency,dim = 2, method="eigen")
      aa <- Mclust(sqrt(n)*Xhat,2,modelNames = c("VVV"))
      error.gmm <- min(sum(abs(aa$classification - tau)),
                            sum(abs(3 - aa$classification - tau)))/n
      
      bb <- kmeans(sqrt(n)*Xhat,2, iter.max = 50)
      error.kmeans <- min(sum(abs(bb$cluster - tau)),
                               sum(abs(3 - bb$cluster - tau)))/n
      
      error.max = max(error.max,error.kmeans)
      
      all <- rbind(all,data.frame(n=c(n,n), err = c(error.gmm,error.kmeans),method=c("gmm","kmeans")))
      
    }
    error.bayes <- bayes.mvn(n, xx[1,], xx[2,], Sigma1, Sigma2, rho[1], rho[2], 1000000)
    error.log <- error.max/(log(nseq[1])/nseq[1])*log(n)/n

    all <- rbind(all, data.frame(n=c(n,n), 
                                 err = c(error.bayes,error.log),
                                 method=c("oracle bayes","log bound")))
  }

  all.stat <- ddply(all,.(n,method),
    function(nn) data.frame(mean=mean(nn$err), se=sd(nn$err)/sqrt(nmc)))
  ggplot(all.stat)+aes(x=n,y=mean,ymin=mean-2*se,ymax=mean+2*se,color=method)+
    geom_line()+
    geom_ribbon(aes(linetype=NA),alpha=.2)+
    ylab("mean error rate")+
    scale_y_log10()
  ggsave("../Figures/gmm_kmeans_bayes.pdf", width=6,height=4)
  
  all
}

experiment3.dan <- function(nmc=100, nseq=seq(from = 1000, to = 4000, by = 250)){
  start.time <- Sys.time()
  rho <- c(0.6,0.4)
  B <- matrix(c(0.42,0.42,0.42,0.5), nrow = 2, ncol = 2)
  x <- eigen(B)
  xx <- x$vectors %*% diag(sqrt(x$values))
  
  all <- data.frame(n=numeric(0),error=numeric(0),name=character(0))
  error.max = 0
  for(n in nseq){
    cat('n=',n,'\n')
    for(j in 1:nmc){
      A <- rg.SBM(n, B, rho,condition=TRUE)
      tau <- A$tau
      X <- xx[tau,]
      
      Xhat <- stfp(A$adjacency,dim = 2, method="eigen")
      
      W <- procrustes(Xhat,X)
      sq.error <- norm(Xhat%*%W - X,"f")^2
      all = rbind(all,data.frame(n=n,sq.error=sq.error,name="Observed"))
      all = rbind(all,data.frame(n=n,sq.error=sqrt(expectation.apv(Xhat)),name="Predicted"))
      cat('\rMC = ',j,' ',format(Sys.time()-start.time))
    }
    all = rbind(all,data.frame(n=n,sq.error=sqrt(expectation.apv(X)),name="Expected"))
    cat('\n')

    p<- ggplot(all)+aes(factor(n),sq.error)+geom_violin()
    ggsave("L2_Error_Violin.pdf",width=6,height=4)
  }

  all
}

expectation.apv <- function(X){
  n <- dim(X)[1]
  SP = diag(t(X) %*% X)
  V = X %*% diag(SP**(-1/2))
  P = X %*% t(X)
  VP =  P*(1-P)-diag(diag(P*(1-P)))
  E <- sum(VP %*% diag(X %*% t(X)))
  for(i in 1:n){
    E <- E+norm(as.matrix(X[i,]))**6 # * norm(as.matrix(V[i,]))**2
#     for(j in 1:n){
#       if(j != i){
#         E <- E+VP[i,j]*norm(as.matrix(X[j,]))**2
#       }
#     }
  }
  E / min(SP**2)
}

asge = function(A,d){
S = eigen(A)
Xhat = S$vectors[,1:d] %*% diag(sqrt(S$values[1:d]))
}      

binomial.test <- function(p,q,r,s,n1,n2,nmc){
    
    X1 <- rbinom(0.6*nmc,n1,p)
    X2 <- rbinom(0.6*nmc,n2,q)

    Y1 <- rbinom(0.4*nmc,n1,r)
    Y2 <- rbinom(0.4*nmc,n2,s)

    Z1 <- c(X1,Y1)
    Z2 <- c(X2,Y2)
    labels <- c(rep(1,0.6*nmc),rep(-1,0.4*nmc))

    tt1 <- Z1*log(p) + (n1-Z1)*log(1-p) + Z2*log(q) + (n2-Z2)*log(1 - q)
    tt2 <- Z1*log(r) + (n1-Z1)*log(1-r) + Z2*log(s) + (n2-Z2)*log(1 - s)

    tt <- tt1 - tt2

    tt[abs(tt) <= 1e-14] <- -1
    yy <- sign(tt)

    error.bayes <- 1 - sum(sign(yy) == labels)/(nmc)
    return(error.bayes)

}


