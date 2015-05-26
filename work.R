
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
##loadData('SSMindependent')  ## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05
loadData('SSMcorrelated')   ## true values: a=.95, b=1, sigPN=0.2, sigOE=0.05
#### SSMcorrelated (i=2) right answers: a=0.922695, b=1.561002, logL=19.09538
#### SSMcorrelated (i=3) right answers: a=0.9250784, b=1.5130347, sigPN=0.1879640, logL=19.4043
#### SSMcorrelated (i=4) right answers: a=?, b=?, sigPN=?, sigPN=?, logL=?
##i <- 2;   MLE <- c(0.922695,  1.561002)
##i <- 3;   MLE <- c(0.9250784, 1.5130347, 0.1879640)
i <- 4;   MLE <- c(0.95, 1.2, 0.2, 0.1)   ### NOTE: green point (MLE value) made up
self <- pfLL(Rmodel, latent, param[1:i], lower[1:i], upper[1:i], init[1:i], trans[1:i], MLE=MLE)
##myKF <- function(a,b) KF_ll(list(y=self$Cmodel$y, sigOE=0.05, sigPN=0.2, a=a, b=b))

self$initialPSOInvestigation()

message('START')
self$iterMVNapproxQuadLMfit()
self$iterMVNapproxQuadLMfit()
self$iterMVNapproxQuadLMfit()
self$iterMVNapproxQuadLMfit()
message('FINISH')








## using ppr()

eigen(self$quadLMABC$A)

## never ut this function into use
fitPPR = function() {
    PPR <- ppr(self$x, matrix(self$y), nterms = self$d, max.terms = self$d)
    print(summary(PPR))
    par(mfrow = c(self$d,1))
    dev.new()
    plot(PPR)
},

d <- 2
V <- array(c(1,-1,0,1), c(2,2))
n <- 100
xlist <- lapply(1:d, function(i) (1-n/2):(n/2))
X <- as.matrix(expand.grid(xlist)) / 100
colnames(X) <- paste0('x', 1:d)
Xp <- X %*% V
##y <- sin(Xp[,1]*pi) + (Xp[,2])^2
y <- -3*(Xp[,1])^3##-3*(Xp[,1]-50)^2 + abs(Xp[,2])

nterms <- 2
##pp <- ppr(cbind(x1,x2), matrix(y), nterms = nterms, max.terms = 5, sm.method='gcvspline')
pp <- ppr(X, matrix(y), nterms = nterms, max.terms = 4)
summary(pp)
par(mfrow = c(nterms,1))
plot(pp)
#plot(pp$gofn, type='b')

## examples from ppr() method:

attach(rock)
area1 <- area/10000; peri1 <- peri/10000
rock.ppr <- ppr(log(perm) ~ area1 + peri1 + shape, data = rock, nterms = 2, max.terms = 5)
rock.ppr
summary(rock.ppr)
par(mfrow = c(3,2))   # maybe: , pty = "s")
plot(rock.ppr, main = "ppr(log(perm)~ ., nterms=2, max.terms=5)")
plot(update(rock.ppr, bass = 5), main = "update(..., bass = 5)")
plot(update(rock.ppr, sm.method = "gcv", gcvpen = 2),
     main = "update(..., sm.method=\"gcv\", gcvpen=2)")
cbind(perm = rock$perm, prediction = round(exp(predict(rock.ppr)), 1))
detach()

## plot the x-points generated from PSO and from MVN
dev.new()
par(mfrow = c(1,1))
plot(self$psoX)
points(self$mvnX, col='red')

## fit LM and plot using MVN points
self$x<-self$mvnX; self$y<-self$mvnY; self$v<-self$mvnV
self$fitModel()
self$plot()

## fit LM and plot using PSO points
self$x<-self$psoX; self$y<-self$psoY; self$v<-self$psoV
self$fitModel()
self$plot()





## finding the difference between logL from KF_ll() and self$CpfLLnf()

myKF <- function(a,b) KF_ll(list(y=data$y, sigOE=0.05, sigPN=0.2, a=a, b=b))
self$CpfLLnf$setM(500000)
self$CpfLLnf$setM(100000)


x1 <- c(0.942619, 1.156373)  ## this is what pfLL() algo finds always (SSMcorrelated, i=2)
x2 <- c(0.8871102, 2.2625421) ## this is the maximum that KF_ll() aka mfKF() finds

thisXpoint <- x2   ## just change this!!

myKF(thisXpoint[1], thisXpoint[2])
self$CpfLLnf$run(thisXpoint)

pfll <- numeric()
it <- 10
for(i in 1:it) {
    pfll[i] <- self$CpfLLnf$run(thisXpoint)
}
pfll
mean(pfll)
var(pfll)
varll <- self$CpfLLnf$getV()
varll
mean(varll)
var(varll)

## correct L(y1), since they are different....
sigPN <- 0.2
sigOE <- 0.05
varPN <- sigPN^2
varOE <- sigOE^2
y <- data$y
a <- thisXpoint[1]
b <- thisXpoint[2]
mu_x <- b / (1-a)
var_x <- varPN / (1-a^2)
mu_y <- mu_x
var_y <- var_x + varOE
l1 <- dnorm(y[1], mu_y, sqrt(var_y), log=TRUE)
l1



## figuring out why the surfaces are fitting poorly
## some really cool 3D plots in here!

rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
model <- 'SSMcorrelated'
loadData(model)
library(plot3D)

myKF <- function(a,b) KF_ll(list(y=data$y, sigOE=0.05, sigPN=0.2, a=a, b=b))

a <- seq(0.86,0.98,by=0.002)
b <- seq(0.4,2.8,by=0.02)
M <- mesh(a, b)
x <- M$x
y <- M$y
ll <- array(NA, dim(x))
for(i in 1:dim(x)[1]) for(j in 1:dim(x)[2]) ll[i,j] <- myKF(a[i],b[j])
##for(i in 1:dim(x)[1]) for(j in 1:dim(x)[2]) ll[i,j] <- self$CpfLLnf$run(c(mu[i],b[j]))
llsave <- ll

minn <- 18
ll[ll<minn] <- minn
surf3D(x, y, ll, phi=0, theta=0, xlab='a', ylab='b')

for(th in seq(0, 180, by=30))
    surf3D(x, y, ll, phi=0, theta=th, main = th, xlab='a', ylab='b')

imax <- which(ll>minn)
cbind(x[imax], y[imax], ll[imax])

maxx <- max(ll)
maxx
minn
trun <- pchisq(2*(maxx-minn), 2)
##12.5% !!!  high-end 30%   middle 20% ?  ## still looks good at 65%
trun

ll <- llsave
##keepInd <- 2*(max(y)-y) < qchisq(self$trunc, self$d)  ## from pfLL.R
discardInd <- which(2*(max(ll)-ll) >= qchisq(trun, 2))
ll[discardInd] <- minn
surf3D(x, y, ll, phi=0, theta=0, xlab='a', ylab='b')

## now try fitting quadratic LM to these heavily reduced points
yLM <- as.numeric(llsave)
n <- length(yLM)
xDF <- expand.grid(a, b)
xLM <- array(NA, dim(xDF))
xLM[1:n, 1] <- xDF[1:n, 1]
xLM[1:n, 2] <- xDF[1:n, 2]
class(xLM)
dimnames(xLM) <- list(NULL, c('a','b'))

keepInd <- yLM>minn
yLM <- yLM[keepInd]
xLM <- xLM[keepInd,]

xterms <- c(c('a','b'), paste0('I(', c('a','b'), '^2)'), if(2 > 1) combn(c('a','b'), 2, function(x) paste0(x, collapse=':')) else character())
form <- as.formula(paste0('y ~ ', paste0(xterms, collapse=' + ')))
df <- as.data.frame(cbind(xLM, y = yLM))
fittedModel <- lm(form, data = df)
summary(fittedModel)
coef <- fittedModel$coef
coef
A <- array(0, c(2, 2), dimnames = list(c('a','b'), c('a','b')))
for(i in seq_along(c('a','b'))) {
    A[i, i] <- coef[paste0('I(', c('a','b')[i], '^2)')]
    if(i!=2) for(j in (i+1):2) A[i, j] <- A[j, i] <- 1/2 * coef[paste0(c('a','b')[i], ':', c('a','b')[j])]
}
b <- array(coef[c('a','b')], c(2, 1), dimnames = list(c('a','b'), NULL))
c <- coef['(Intercept)']
if(any(A==0) || any(b==0) || (c==0) ) stop('coefficients not extracted correctly')
if(any(eigen(A)$values > 0)) warning('Negative-quadratic model fit is NOT strictly concave down')
eigen(A)
paramT <- t(-1/2 * solve(A) %*% b)[1, ]
logL <- (-1/4 * t(b) %*% solve(A) %*% b + c)[1,1]
paramT
##        a         b 
##0.9209136 1.5967389 
logL
##[1] 19.09287

## let's see how these LM predictions compare to our actual data
newMin <- 18.9
imax <- which(ll>newMin)
length(imax)
xTop <- x[imax]
yTop <- y[imax]
llTop <- ll[imax]
ix <- sort(llTop, index.return=TRUE)$ix
xTop <- xTop[ix]
yTop <- yTop[ix]
llTop <- llTop[ix]
ddd <- cbind(xTop, yTop, llTop)
ddd
apply(ddd, 2, mean)
##    xTop     yTop    llTop 
## 0.91900  1.64000 18.97373 
## answer: very well!

optimInit <- c(0.942619, 1.156373)
myKFoptim <- function(ab) myKF(ab[1], ab[2])
optimOut <- optim(optimInit, myKFoptim, control=list(fnscale=-1))
optimOut$convergence
optimOut$par
##[1] 0.922695 1.561002
optimOut$value
##[1] 19.09538

## here's for SSMcorrelated (i=3) with sigPN on log-scale
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
model <- 'SSMcorrelated'
loadData(model)
library(plot3D)
myKF3 <- function(a,b,logsigPN) KF_ll(list(y=data$y, sigOE=0.05, sigPN=exp(logsigPN), a=a, b=b))
optimInit3 <- c(0.922695, 1.561002, log(0.2))
myKFoptim3 <- function(abLPN) myKF3(abLPN[1], abLPN[2], abLPN[3])
optimOut3 <- optim(optimInit3, myKFoptim3, control=list(fnscale=-1))
optimOut3$convergence
par3 <- optimOut3$par
c(par3[1], par3[2], exp(par3[3]))
## 0.9250784 1.5130347 0.1879640
optimOut3$value
##[1] 19.4043

## developing MVN simulation for a new cloud of points


self$bestX
apply(self$x, 2, mean)
cov(self$x)

n <- 400
xNew <- rmvnorm(n, self$bestX, cov(self$x))
for(i in 1:n) self$CpfLLnf$run(xNew[i, ])

par(mfrow = c(1,1))
plot(self$x)
points(xNew, col='red')



## comparing pfLL performance to KF, at true values, maximized value, etc..

self$CpfLLnf
KF_ll(c(as.list(self$max$param), list(y=data$y)))

KF_ll(c(data,inits))  ## [1] 19.33315
self$param
self$CpfLLnf$run(c(20, 1, log(0.2), log(0.05)))
## MLE values (according to optim()): mu=20, b=2.34, sigPN=0.193, sigOE=0.0001
## MLE log-likelihood (according to optim()): 21.75809



## analytically finding neg-quadratic surface peak,
## *and* check that fitted surface is strictly concave down.

coef <- self$fittedModel$coef
A <- array(0, c(self$d, self$d), dimnames = list(self$param, self$param))
for(i in seq_along(self$param)) {
    A[i, i] <- coef[paste0('I(', self$param[i], '^2)')]
    if(i!=self$d) for(j in (i+1):self$d) A[i, j] <- A[j, i] <- 1/2 * coef[paste0(self$param[i], ':', self$param[j])]
}
b <- array(coef[self$param], c(self$d, 1), dimnames = list(self$param, NULL))
c <- coef['(Intercept)']

self$param
coef
A
b
c

eigen(A)$values   ## should be all negative ???

-1/2 * solve(A) %*% b
t(t(self$max$paramT))

-1/4 * t(b) %*% solve(A) %*% b + c
self$max$logL



## testing accuracy of PF LL variance estimate

rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMindependent')

Rpf <- buildPF(Rmodel, 'x')
Cmodel <- compileNimble(Rmodel)
Cpf <- compileNimble(Rpf, project = Rmodel)

n <- 500
ll <- numeric()
varll <- numeric()

for(i in 1:n) {
    print(i/n)
    ll[i] <- Cpf$run(500)
    varll[i] <- Cpf$getVarLL()
}

ll
mean(ll)
var(ll)

varll
mean(varll)
var(varll)



## testing how to populate & resize NIMBLE vectors / arrays

library(nimble)

nfDef <- nimbleFunction(
    setup = function() {
        d <- 4
        increment <- 10  ## must be > 1
        x <- array(as.numeric(NA), c(increment, d))
        y <- rep(as.numeric(NA), increment)
        cur <- 0
        max <- increment
    },
    run = function() {
        if(cur == max) {
            max <<- max + increment
            print('max increased to: ', max)
            xTemp <- x
            yTemp <- y
            setSize(x, max, d)
            setSize(y, max)
            x[1:cur, 1:d] <<- xTemp
            y[1:cur]      <<- yTemp
        }
        cur <<- cur + 1
        print('cur increased to: ', cur)
        for(i in 1:d) x[cur, i] <<- 10*cur + i
        y[cur] <<- cur/10
    },
    methods = list(
        getX = function() {
            returnType(double(2))
            return(x[1:cur, 1:d])
        },
        getY = function() {
            returnType(double(1))
            return(y[1:cur])
        }
    )
)

Rnf <- nfDef()
Cnf <- compileNimble(Rnf)

up <- 35
for(i in 1:up) Cnf$run()

Cnf$getX()
Cnf$getY()


