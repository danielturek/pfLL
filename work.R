

rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMindependent')

paramNodes <- c('mu')
lower <- c(19)
upper <- c(22)

paramNodes <- c('mu', 'b')
lower <- c(19, 0)
upper <- c(22, 10)

pfLLobj <- pfLLClass$new(Rmodel, latent, paramNodes, 1000)
set.seed(0);     pfLLobj$run(c(20.5, 1)) # 20.32805

control <- list(fnscale=-1, trace=1, trace.stats=TRUE, REPORT=1, vectorize=TRUE, maxit=20)
set.seed(0); system.time(out <- psoptim(paramNodes, pfLLobj$run, lower=lower, upper=upper, control=control))
out$par    # 20.216061,  1.450972
out$value  # -22.38359

trace <- extractPSOtrace(out)

plotPSOtrace(trace)







## running of a PF on param values
## x <- seq(20, 20.5, by=0.05)
## ll <- sapply(x, CpfLL$run)
## plot(x, ll)
## max(ll) ## 20.45735
## x[which(ll==max(ll))] ## 20.15


## build & run PF
## Rpf <- buildPF(Rmodel, latent)
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


