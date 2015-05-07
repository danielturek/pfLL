
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMindependent')
## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05


paramNodes <- c('mu'); lower <- c(19); upper <- c(22)
## psoIt = 3
## [1] 20.11637
## [1] -21.22038
## Multiple R-squared:  0.9222,	Adjusted R-squared:  0.9175

paramNodes <- c('mu','b'); lower <- c(19,0); upper <- c(22,10)
## psoIt = 3
## [1] 20.142429  2.255136
## [1] -21.04907
## Multiple R-squared:  0.7848,	Adjusted R-squared:  0.7489

paramNodes <- c('mu','b','sigPN','sigOE'); lower <- c(19,0,.01,.01); upper <- c(22,4,1,1)
## psoIt = 3
## [1] 19.1654744  0.9124417  0.2935520  0.1689781
## [1] 12.74345
## Multiple R-squared:  0.6636,	Adjusted R-squared:  0.4825
## psoIt = 50:
## [1] 20.08460637  1.84003807  0.18720946  0.02067147
## [1] -23.55684
## Multiple R-squared:  0.1718,	Adjusted R-squared:  0.1547


self <- pfLL(Rmodel, latent, paramNodes, lower, upper, m=5000, psoIt=50)
print(self$psoOut$par); print(self$psoOut$value)



## now write some spline fitting!!!!

x <- self$psoTrace$x
y <- self$psoTrace$y
n <- self$psoTrace$n
p <- self$psoTrace$p
paramNodes
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


