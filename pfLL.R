
library(nimble)
library(R6)
model <- 'SSMindependent'
load(paste0('~/GitHub/autoBlock/data/model_', model, '.RData'))
Rmodel <- nimbleModel(code, constants, data, inits)
latent <- 'x'

pfLL <- function(model, latent, param, m = 10000, rep = 1, makePlot = TRUE) {
    pfLLobj <- pfLLClass$new()
    pfLLobj$initModel(model)
    pfLLobj$initPF(latent)
    pfLLobj$initParam(param)
    pfLLobj$initRep(rep)
    pfLLobj$calcLL(m)
    if(makePlot) pfLLobj$makePlot()
    pfLLobj
}; pfLLClass <- R6Class(
    'pfLLClass',
    public = list(
        Rmodel = NULL,
        Cmodel = NULL,
        initModel = function(modelarg) {
            md <- modelarg$modelDef
            self$Rmodel <- md$newModel()
            nimCopy(from = modelarg, to = self$Rmodel, logProb = TRUE)
            for(var in ls(modelarg$isDataEnv)) self$Rmodel$isDataEnv[[var]] <- modelarg$isDataEnv[[var]]
            self$Cmodel <- compileNimble(self$Rmodel)
        },
        Rpf = NULL,
        Cpf = NULL,
        initPF = function(latent) {
            self$Rpf <- buildPF(self$Rmodel, latent)
            self$Cpf <- compileNimble(self$Rpf, project = self$Rmodel)
        },
        param = NULL,
        paramNames = NULL,
        nParamValues = NULL,
        initParam = function(param) {
            self$param <- param
            self$paramNames <- dimnames(param)[[2]]
            self$nParamValues <- dim(param)[1]
            if(!all(self$paramNames %in% self$Cmodel$getNodeNames(returnScalarComponents = TRUE)))
                stop('problem with param names')
        },
        rep = NULL,
        initRep = function(rep) {
            self$rep <- rep
        },
        ll = NULL,
        calcLL = function(m) {
            self$ll <- array(NA, c(self$nParamValues, self$rep))
            for(i in 1:self$nParamValues) {
                cat(paste0(round((i-1)/self$nParamValues*100,0), '%\n'))
                self$setModelParamValues(row = i)
                for(j in 1:self$rep)     self$ll[i, j] <- self$Cpf$run(m)
            }
        },
        setModelParamValues = function(row) {
            values(self$Cmodel, self$paramNames) <- as.numeric(self$param[row, ])
        },
        makePlot = function() {
            switch(length(self$paramNames),
                   plot(rep(self$param[,1],self$rep),
                        as.numeric(self$ll),
                        xlab=self$paramNames,
                        ylab='log-likelihood estimate'),  ## plot for 1 param
                   { require(plot3D)
                     scatter3D(rep(self$param[,1],self$rep),
                               rep(self$param[,2],self$rep),
                               as.numeric(self$ll),
                               xlab=self$paramNames[1],
                               ylab=self$paramNames[2],
                               zlab='log-likelihood estimate') }  ## plot for 2 params
                   )
        }
    )
)


## 1 param
param <- expand.grid(list(mu = seq(19, 22, by=0.2)))
out <- pfLL(Rmodel, latent, param, rep=5, m=3000)

## 2 params
param <- expand.grid(list(mu = seq(19.6, 21, by=0.2), b = seq(-2, 5, by=1)))
out <- pfLL(Rmodel, latent, param, rep=3, m=2000)








