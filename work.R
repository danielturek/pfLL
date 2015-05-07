
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMindependent')
## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05


paramNodes <- c('mu'); lower <- c(19); upper <- c(22)

paramNodes <- c('mu','b'); lower <- c(19,-10); upper <- c(22,10)

paramNodes <- c('mu','b','sigPN','sigOE'); lower <- c(19,-10,.01,.01); upper <- c(22,4,1,1)
## psoIt = 20:
## [1] 20.18847064  1.56384181  0.20058355  0.02346878
## [1] -22.13378

self <- pfLL(Rmodel, latent, paramNodes, lower, upper, m=5000, psoIt=20)
print(self$psoOut$par); print(self$psoOut$value)





yPred <- as.numeric(predict(self$quadLM))



iParam <- 4

plot(self$x[,iParam], self$y, xlab=self$paramNodes[iParam], ylab='log-likelihood')
sortedX <- sort(self$x[,iParam], index.return = TRUE)
thisX <- sortedX$x
thisY <- yPred[sortedX$ix]
lines(thisX, thisY, col='blue')
hullInd <- chull(thisX, thisY)
while(hullInd[1]!=min(hullInd)) hullInd <- c(hullInd[2:length(hullInd)], hullInd[1])
hullInd <- hullInd[1:which(hullInd==max(hullInd))]
hullX <- thisX[hullInd]
hullY <- thisY[hullInd]
lines(hullX, hullY, col='red', lwd=2)



## now write some spline fitting!!!!
require(mgcv)








## running of a PF on param values
## x <- seq(20, 20.5, by=0.05)
## ll <- sapply(x, CpfLL$run)
## plot(x, ll)
## max(ll) ## 20.45735
## x[which(ll==max(ll))] ## 20.15


## build & run PF
## Rpf <- buildPF(Rmodel, latent)
## Cmodel <- compileNimble(Rmodel)
## Cpf <- compileNimble(Rpf, project = Rmodel)
## set.seed(0)
## Cpf$run(10000)


## Creating data file of LLs, 2D and 3D plots
## ## 1 param
## param <- expand.grid(list(mu = seq(19, 22, by=0.2)))
## out1 <- pfLL(Rmodel, latent, param, rep=5, m=3000)
## ## 2 params
## param <- expand.grid(list(mu = seq(19.6, 21, by=0.2), b = seq(-5, 7, by=2)))
## out2 <- pfLL(Rmodel, latent, param, rep=3, m=2000)
## ## save for spline testing!
## save(out1, out2, file='~/GitHub/pfLL/pfLL.RData')
##
## rm(list=ls())
## load('~/Github/pfLL/pfLL.RData')
## require(plot3D)
## plot(out1$x[, 1], out1$y, xlab=names(out1$param), ylab='ll')  ## 1 param plot
## scatter3D(out2$x[, 1], out2$x[, 2], out2$y, xlab=names(out2$param)[1], ylab=names(out2$param)[2], zlab='ll')  ## 2 param plot
## max(out2$y) ## 20.51823
## out2$x[which(out2$ll==max(out2$ll)), ] ## mu=20, b=3


