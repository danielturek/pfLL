
library(nimble)
require(R6)
model <- 'SSMindependent'
load(paste0('~/GitHub/autoBlock/data/model_', model, '.RData'))
Rmodel <- nimbleModel(code, constants, data, inits)
latent <- 'x'

pfLL <- function(model, latent, param, m = 10000, rep = 1, makePlot = TRUE) {
    require(R6)
    require(nimble)
    pfLLobj <- pfLLClass$new()
    pfLLobj$initModel(model)
    pfLLobj$initPF(latent)
    pfLLobj$initParam(param)
    pfLLobj$calcLL(m, rep)
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
        ll = NULL,
        calcLL = function(m, rep) {
            self$ll <- array(NA, c(self$nParamValues, rep))
            for(i in 1:self$nParamValues) {
                self$setModelParamValues(row = i)
                print(self$Cmodel[[self$paramNames[1]]])    ## temporary testing
                for(j in 1:rep)     self$ll[i, j] <- self$Cpf$run(m)
            }
        },
        setModelParamValues = function(row) {
            values(self$Cmodel, self$paramNames) <- as.numeric(self$param[row, ])
        },
        makePlot = function() {
            switch(length(self$paramNames),
                   plot(x=rep(self$param[,1], dim(self$ll)[2]), y = self$ll, xlab=self$paramNames, ylab='log-likelihood estimate'),  ## plot for 1 param
                   print('need to implement surface plotting')  ## plot for 2 params
                   )
        }
    )
)

param <- data.frame(mu = seq(19, 21, by=0.2))

pfLL(Rmodel, latent, param, rep=5)








