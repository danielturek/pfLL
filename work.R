
##rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
##loadData('SSMindependent')  ## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05
loadData('SSMcorrelated')   ## true values: a=.95, b=1, sigPN=0.2, sigOE=0.05
#### SSMcorrelated (i=2) right answers: a=0.922695, b=1.561002, logL=19.09538
#### SSMcorrelated (i=3) right answers: a=0.9250784, b=1.5130347, sigPN=0.1879640, logL=19.4043
#### SSMcorrelated (i=4): a=0.9191, b=1.6333, sigPN=0.1955, sigOE=1.476e-5, logL=20.38073
i <- 2
self <- pfLL(Rmodel, latent, param[1:i], lower[1:i], upper[1:i], trans[1:i], trunc=0.6, run=FALSE)
cnf <- self$CpfLLnf

## seeing what Nelder Mead optimization does with this

myKF <- function(a,b) KF_ll(list(y=data$y, sigOE=0.05, sigPN=0.2, a=a, b=b))
NM <- function(...) NMClass$new(...)
NMClass <- R6Class(
    'NMClass',
    public = list(
        cnf = NULL,
        x = NULL, y = NULL, v = NULL, n = NULL,
        xhist = NULL,
        alpha = NULL, beta = NULL, gamma = NULL, delta = NULL, k = NULL, KF = NULL,
        initialize = function(xinit, cnf, k = 1.25, m = 5000, KF = FALSE) {
            set.seed(0)
            self$cnf <- cnf
            self$k <- k
            self$KF <- KF
            self$cnf$setM(m)
            cnf$reset()
            optim(xinit, cnf$run, control = list(fnscale = -1, maxit = 13))
            ii <- c(1, 5, 6)
            self$x <- cnf$getX()[ii,]
            self$y <- cnf$getY()[ii]
            self$v <- cnf$getV()[ii]
            self$n <- length(self$v)
            self$sortXYV()
            self$xhist <- self$x
            self$alpha <- 1
            self$beta <- 0.5
            self$gamma <- 2
            self$delta <- 0.9##0.5
        },
        sortXYV = function() {
            ix <- sort(self$y, decreasing=TRUE, index.return=TRUE)$ix
            self$x <- self$x[ix,]
            self$y <- self$y[ix]
            self$v <- self$v[ix]
        },
        incPrec = function() {
            message('increasing precision')
            cnf$scaleM(self$k)
            message('m is now: ', cnf$getM())
        },
        decPrec = function() {
            message('decreasing precision')
            cnf$scaleM(1/self$k)
            message('m is now: ', cnf$getM())
        },
        calcOLD = function(xcalc) {
            x <- cnf$getX()
            y <- cnf$getY()
            v <- cnf$getV()
            for(i in 1:dim(x)[1]) {
                if(all(xcalc == x[i,])) {
                    self$xhist <- rbind(self$xhist, xcalc)
                    return(c(y[i], v[i]))
                }
            }
            test <- cnf$run(xcalc)
            if(test == -Inf) return(c(-Inf, Inf))
            print(xcalc)
            print(test)
            stop('coundlt find right xcalc point...')
        },
        calc = function(xcalc) {
            if(self$KF) {
                ll <- myKF(xcalc[1], xcalc[2])
                return(c(ll, 1))
            }
            ll <- cnf$run(xcalc)
            if(ll == -Inf) return(c(-Inf, Inf))
            lastY <- cnf$getLastY(0)
            lastV <- cnf$getLastV(0)
            if(ll != lastY) stop('wrong')
            return(c(lastY, lastV))
        },
        accept = function(xacc, yacc, vacc, pos) {
            self$x[pos,] <- xacc
            self$y[pos] <- yacc
            self$v[pos] <- vacc
            self$xhist <- rbind(self$xhist, xacc)
            self$sortXYV()
        },
        adapt = function() {
            T <- sum(((self$y - mean(self$y)) / sqrt(self$v))^2)
            if(!self$KF) message('T is: ', T) else message('NA')
            pn <- pchisq(T, df=self$n)
            if(!self$KF) message('pn is: ', pn)  else message('NA')
            if(pn < 0.5) {  ## fail to rejust H0: means are equal
                self$incPrec()
                ##if(!self$KF) message('fail to reject H0; would increase precision here') else message('NA')
            } else {
                self$decPrec()
                ##if(!self$KF) message('reject H0; decrease precision here') else message('NA')
            }
        },
        step = function() {
            self$adapt()
            for(i in 1:1) {
                thisx <- self$x[i,]
                thisfv <- self$calc(thisx)
                self$accept(thisx, thisfv[1], thisfv[2], i)
                message('recalc xbest: ', thisfv[1])
                print(thisx)
            }
            self$sortXYV()
            xbest <- self$x[1,]
            xsecond <- self$x[2,]
            y1 <- self$y[1]
            y2 <- self$y[2]
            message('xbest: ', y1)
            print(xbest)
            message('xsecond: ', y2)
            print(xsecond)
            message('xworst: ', self$y[3])
            print(self$x[3,])
            xc <- apply(self$x[-self$n,], 2, mean)
            message('xc:')
            print(xc)
            xworst <- self$x[self$n,]
            yworst <- self$y[self$n]
            xr <- xc + self$alpha*(xc-xworst)
            ##messge('calculate xr')
            fvr <- self$calc(xr)
            fr <- fvr[1]
            message('calc xr: ', fr)
            print(xr)
            if(y1 >= fr & fr > y2) { message('reflection'); self$accept(xr, fr, fvr[2], self$n); return() }
            if(fr > y1) {
                xe <- xc + self$gamma*(xr-xc)
                fve <- self$calc(xe)
                fe <- fve[1]
                message('calc xe: ', fe)
                message('variance of xe calculation: ', fve[2])
                print(xe)
                ##if(fe > fr)  { self$accept(xe, fe, fve[2], self$n); return() }  ## greedy min.
                if(fe > y1)  { message('expansion'); self$accept(xe, fe, fve[2], self$n); return() }  ## greedy exp.
                message('reflection second chance'); self$accept(xr, fr, fvr[2], self$n); return()
            }
            if(fr <= y2) {
                if(fr > yworst) {
                    xcon <- xc + self$beta * (xr - xc)
                    fvcon <- self$calc(xcon)
                    fcon <- fvcon[1]
                    message('calc xcontract 2: ', fcon)
                    print(xcon)
                    ##if(fcon >= fr) { message('contraction 2 from xr'); self$accept(xcon, fcon, fvcon[2], self$n); return() }
                    ## DO NOT ALLOW SHRINKS
                    if(fcon < fr) { message('would have shrunk, not going to') }
                    if(TRUE) { message('contraction 2 from xr'); self$accept(xcon, fcon, fvcon[2], self$n); return() }
                }
                if(fr <= yworst) {
                    xcon <- xc + self$beta * (xworst - xc)
                    fvcon <- self$calc(xcon)
                    fcon <- fvcon[1]
                    message('calc xcontract 1: ', fcon)
                    print(xcon)
                    ##if(fcon > yworst) { message('contraction 1 from xworst'); self$accept(xcon, fcon, fvcon[2], self$n); return() }
                    ## DO NOT ALLOW SHRINKS
                    if(fcon <= yworst) { message('would have shrunk, not going to') }
                    if(TRUE) { message('contraction 1 from xworst'); self$accept(xcon, fcon, fvcon[2], self$n); return() }
                }
            }
            message('shrink')
            for(j in 2:self$n) {
                xnew <- xbest + self$delta * (self$x[j,] - xbest)
                fvnew <- self$calc(xnew)
                fnew <- fvnew[1]
                self$accept(xnew, fnew, fvnew[2], j)
            }
        }
    )
)
xinit <- c(0.866, 2.70)

myNM <- NM(xinit, cnf, KF=FALSE, k=1.1)
myNMKF <- NM(xinit, cnf, KF=TRUE)
iter <- 20
myNM$cnf$getM()
myNM$cnf$setM(50000)
for(i in 1:iter) { message(i); myNM$step() }
for(i in 1:iter) { message(i); myNMKF$step() }

x <- myNM$xhist
ind <- c(30:40)
#ind <- 1:dim(x)[1]
rng <- apply(x[ind,], 2, range)
a <- seq(rng[1,1],rng[2,1],length.out=100)
b <- seq(rng[1,2],rng[2,2],length.out=100)
myKF <- function(a,b) KF_ll(list(y=data$y, sigOE=0.05, sigPN=0.2, a=a, b=b))
ll <- outer(a, b, FUN = myKF)
plotpnts <- function(ind, x, xoffset=0, ...) text(x[ind,1]+xoffset, x[ind,2], labels = ind, cex=0.5, ...)
contour(a, b, ll, nlevels=20)
points(0.922695, 1.561002, pch=19, col='green')  ## plots the MLE
plotpnts(ind, myNM$xhist,    col='red', xoffset=0.00002)
##plotpnts(ind, myNMKF$xhist, col='green')

myNM$x
myNM$y
myNM$v

cnf$run(c(0.849125, 2.725312))





## examining plots of psoXYV and mvnXYV, and fitting models

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

##rm(list=ls())
psource('~/GitHub/pfLL/pfLL.R')
model <- 'SSMcorrelated'
loadData(model)
library(plot3D)

myKF <- function(a,b) KF_ll(list(y=data$y, sigOE=0.05, sigPN=0.2, a=a, b=b))

a <- seq(0.86,0.98,by=0.002)
b <- seq(0.4,2.8,by=0.02)
ll <- outer(a, b, FUN = myKF)
llsave <- ll
M <- mesh(a, b)
x <- M$x
y <- M$y

contour(a, b, ll, nlevels=200)


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
newMin <- 17
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



## finding the "right" answers for SSMcorrelated, using optim()

## SSMcorrelated (i=2)
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMcorrelated')
myKF2 <- function(a,b) KF_ll(list(y=data$y, sigOE=0.05, sigPN=0.2, a=a, b=b))
optimInit2 <- c(0.942619, 1.156373)
myKFoptim2 <- function(ab) myKF2(ab[1], ab[2])
optimOut2 <- optim(optimInit2, myKFoptim2, control=list(fnscale=-1))
optimOut2$convergence
optimOut2$par
##[1] 0.922695 1.561002
optimOut2$value
##[1] 19.09538

## SSMcorrelated (i=3) with sigPN on log-scale
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMcorrelated')
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

## SSMcorrelated (i=4) with sigPN and sigOE on log-scale
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
loadData('SSMcorrelated')
myKF4 <- function(a,b,logsigPN,logsigOE) KF_ll(list(y=data$y, sigOE=exp(logsigOE), sigPN=exp(logsigPN), a=a, b=b))
optimInit4 <- c(0.922695, 1.561002, log(0.2), log(0.05))
myKFoptim4 <- function(abLPNLOE) myKF4(abLPNLOE[1], abLPNLOE[2], abLPNLOE[3], abLPNLOE[4])
optimOut4 <- optim(optimInit4, myKFoptim4, control=list(fnscale=-1))
optimOut4$convergence
par4 <- optimOut4$par
c(par4[1], par4[2], exp(par4[3]), exp(par4[4]))
## 9.190980e-01 1.633285e+00 1.955246e-01 1.475635e-05
optimOut4$value
## [1] 20.38073



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


