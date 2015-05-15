

################
### SSMindependent
################

## better parameterization: mean and autocorrelation
rm(list = ls())
modelName <- 'SSMindependent'
modelfileName <- paste0('~/GitHub/pfLL/data/model_', modelName, '.RData')
library(nimble)
code <- quote({
    mu ~ dnorm(0, sd = 1000)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
    x[1] ~ dnorm(mu, sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    a <- 1-(b/mu)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})
t <- 100
constants <- list(t = t)
Rmodel <- nimbleModel(code = code, constants = constants)
Rmodel$mu <- 1/(1-.95)
Rmodel$b <- 1
Rmodel$sigPN <- .2
Rmodel$sigOE <- .05
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('mu','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data <- list(y = as.numeric(Rmodel$y))
inits <- list(mu = Rmodel$mu, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = as.numeric(Rmodel$x))
latent <- 'x'
param <- c('mu',     'b',      'sigPN',  'sigOE')
lower <- c(19,      -10,       .000001,  .000001)
upper <- c(22,       10,       1,        1)
trans <- c(identity, identity, log,      log)
nParam <- length(param)
save(code, constants, data, inits, latent, param, lower, upper, trans, nParam, file=modelfileName)



################
### SSMcorrelated
################

## parameterization in terms of slope (autocorelation) and intercept,
## which are highly correlated in mixing
rm(list = ls())
modelName <- 'SSMcorrelated'
modelfileName <- paste0('~/GitHub/pfLL/data/model_', modelName, '.RData')
library(nimble)
code <- quote({
    a ~ dunif(-0.9999, 0.9999)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
    x[1] ~ dnorm(b/(1-a), sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})
t <- 100
constants <- list(t = t)
Rmodel <- nimbleModel(code = code, constants = constants)
Rmodel$a <- .95
Rmodel$b <- 1
Rmodel$sigPN <- .2
Rmodel$sigOE <- .05
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('a','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data <- list(y = as.numeric(Rmodel$y))
inits <- list(a = Rmodel$a, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = as.numeric(Rmodel$x))
latent <- 'x'
param <- c('a',      'b',      'sigPN', 'sigOE')
lower <- c(0,        0.5,      .0001,   .0001)
upper <- c(0.99,     1.5,      1,       1)
trans <- c(identity, identity, log,     log)
nParam <- length(param)
save(code, constants, data, inits, latent, param, lower, upper, trans, nParam, file=modelfileName)



