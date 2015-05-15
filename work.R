
rm(list=ls())
source('~/GitHub/pfLL/pfLL.R')
model <- 'SSMindependent'  ## true values: mu=20, b=1, sigPN=0.2, sigOE=0.05
##model <- 'SSMcorrelated'  ## true values: a=.95, b=1, sigPN=0.2, sigOE=0.05
loadData(model)

## MLE values (according to optim()): mu=20, b=2.34, sigPN=0.193, sigOE=0.0001
## MLE log-likelihood (according to optim()): 21.75809
## d <- list(y=data$y, mu=inits$mu, a=1-inits$b/inits$mu, b=inits$b, sigOE=inits$sigOE, sigPN=inits$sigPN, t=constants$t)
## KF_ll(d)  ## [1] 19.33315, only runs for 'SSMindependent' model

i <- 2
self <- pfLL(Rmodel,latent,param[1:i],lower[1:i],upper[1:i],trans[1:i],m=1000, psoIt=10)






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
        p <- 4
        increment <- 10  ## must be > 1
        x <- array(as.numeric(NA), c(increment, p))
        y <- rep(as.numeric(NA), increment)
        cur <- 0
        max <- increment
    },
    run = function() {
        cur <<- cur + 1
        if(cur > max) {
            print('ERROR: bugs in setSize(), so this will fail!')
            print('ERROR: see Github issues #57, #58')
            max <<- max + increment
            setSize(x, max, p)
            setSize(y, max)
        }
        for(i in 1:p) x[cur, i] <<- 10*cur + i
        y[cur] <<- cur/10
    },
    methods = list(
        getX = function() {
            returnType(double(2))
            return(x[1:cur, 1:p])
        },
        getY = function() {
            returnType(double(1))
            return(y[1:cur])
        }
    )
)

Rnf <- nfDef()
Cnf <- compileNimble(Rnf)

up <- 71
for(i in 1:up) Cnf$run()

Cnf$getX()
Cnf$getY()



























