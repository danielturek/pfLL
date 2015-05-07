
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMindependent')  ## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05
d <- list(y=as.numeric(data$y), mu=inits$mu, a=1-inits$b/inits$mu, b=inits$b, sigOE=inits$sigOE, sigPN=inits$sigPN, t=constants$t)

KF_ll(d)  ## [1] 19.33315

m <- 5000
psoIt <- 10
paramNodes <- c('mu'); lower <- c(19); upper <- c(22); trans <- c(identity)
self1 <- self <- pfLL(Rmodel, latent, paramNodes, lower, upper, trans, m=m, psoIt=psoIt)
paramNodes <- c('mu','b'); lower <- c(19,-10); upper <- c(22,10); trans <- c(identity, identity)
self2 <- self <- pfLL(Rmodel, latent, paramNodes, lower, upper, trans, m=m, psoIt=psoIt)
paramNodes <- c('mu','b','sigPN','sigOE'); lower <- c(19,-10,.01,.01); upper <- c(22,4,1,1); trans <- c(identity, identity, log, log)
self4 <- self <- pfLL(Rmodel, latent, paramNodes, lower, upper, trans, m=m, psoIt=psoIt)


Rmodel <- NULL
save(list=ls(), file = '~/GitHub/pfLL/pfLL_SSMindependent.RData')

load('~/GitHub/pfLL/pfLL_SSMindependent.RData')






## write some spline fitting (???)
## require(mgcv)


rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMcorrelated')  ## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05
d <- list(y=as.numeric(data$y), mu=inits$b/(1-inits$a), a=inits$a, b=inits$b, sigOE=inits$sigOE, sigPN=inits$sigPN, t=constants$t)

KF_ll(d)  ## [1] 19.33315

m <- 5000
psoIt <- 10
paramNodes <- c('a'); lower <- c(0.945); upper <- c(0.955); trans <- c(identity)
self1 <- self <- pfLL(Rmodel, latent, paramNodes, lower, upper, trans, m=m, psoIt=psoIt)
paramNodes <- c('a','b'); lower <- c(0.945,-10); upper <- c(0.955,10); trans <- c(identity, identity)
self2 <- self <- pfLL(Rmodel, latent, paramNodes, lower, upper, trans, m=m, psoIt=psoIt)
paramNodes <- c('a','b','sigPN','sigOE'); lower <- c(0.945,-10,.01,.01); upper <- c(0.955,4,1,1); trans <- c(identity, identity, log, log)
self4 <- self <- pfLL(Rmodel, latent, paramNodes, lower, upper, trans, m=m, psoIt=psoIt)


Rmodel <- NULL
save(list=ls(), file = '~/GitHub/pfLL/pfLL_SSMcorrelated.RData')





