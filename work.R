
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
#model <- 'SSMindependent'  ## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05
model <- 'SSMcorrelated'  ## true values: a=.95, b=1, sigPN=0.2, sigOE=0.05
loadData(model)
i <- 2
self <- pfLL(Rmodel, latent, param[1:i], lower[1:i], upper[1:i], trans[1:i])




## trying MVN simulation to create a good spread of points

par(mfrow = c(1,2))
plot(self$x)

self$bestX
apply(self$x, 2, mean)

cov(self$x)

library(mvtnorm)
n <- 1000
xNew <- rmvnorm(n, self$bestX, cov(self$x))
plot(xNew)


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


