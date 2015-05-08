
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMindependent')
## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05
## MLE values (according to optim()): mu=20, b=2.34, sigPN=0.193, sigOE=0.0001
## MLE log-likelihood (according to optim()): 21.75809

d <- list(y=as.numeric(data$y), mu=inits$mu, a=1-inits$b/inits$mu, b=inits$b, sigOE=inits$sigOE, sigPN=inits$sigPN, t=constants$t)
KF_ll(d)  ## [1] 19.33315

i <- 4

paramNodes <- c('mu',     'b',      'sigPN',  'sigOE')
lower <-      c(19,      -10,       .000001,  .000001)
upper <-      c(22,       10,       1,        1)
trans <-      c(identity, identity, log,      log)
s <- list()
for(i in 1:4) {
    pN<-paramNodes[1:i];  low<-lower[1:i];  upp<-upper[1:i];  tr<-trans[1:i]
    s[[i]] <- self <- pfLL(Rmodel, latent, pN, low, upp, tr, m=5000, psoIt=40)
}

Rmodel <- NULL
save(list=ls(), file = '~/GitHub/pfLL/pfLL_SSMindependent.RData')

load('~/GitHub/pfLL/pfLL_SSMindependent.RData')




#######################
## correlated SSM model
#######################

rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMcorrelated')  ## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05

i <- 2

paramNodes <- c('a',      'b',      'sigPN', 'sigOE')
lower <-      c(0,     0.5,      .0001,   .0001)
upper <-      c(0.99,    1.5,      1,       1)
trans <-      c(identity, identity, log,     log)
s <- list()
for(i in 1:4) {
    pN<-paramNodes[1:i];  low<-lower[1:i];  upp<-upper[1:i];  tr<-trans[1:i]
    s[[i]] <- self <- pfLL(Rmodel, latent, pN, low, upp, tr, m=5000, psoIt=40)
}

Rmodel <- NULL
save(list=ls(), file = '~/GitHub/pfLL/pfLL_SSMcorrelated.RData')






summary(self$quadLM)
self$optimOut

df <- data.frame(as.list(self$optimOut$par))
names(df) <- self$paramNodes
df
predict(self$quadLM, df)
m <- self$quadLM
summary(m)
coef(m)

a <- self$optimOut$par[1]
b <- self$optimOut$par[2]

X <- c(1, a, b, a^2, b^2, a*b)
X * coef(m)

require(plot3D)
x <- self$x[, 1]
y <- self$x[, 2]
z <- self$y
scatter3D(x, y, z, theta = 130)










