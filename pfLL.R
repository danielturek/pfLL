
require(nimble)
require(R6)
require(pso)

loadData <- function(model) {
    load(paste0('~/GitHub/pfLL/data/model_', model, '.RData'), envir=parent.frame())
    eval(quote(Rmodel <- nimbleModel(code, constants, data, inits)), envir = parent.frame())
}

replicateModel <- function(modelarg) {
    md <- modelarg$modelDef
    Rmodel <- md$newModel()
    nimCopy(from = modelarg, to = Rmodel, logProb = TRUE)
    for(var in ls(modelarg$isDataEnv)) Rmodel$isDataEnv[[var]] <- modelarg$isDataEnv[[var]]
    return(Rmodel)
}

## a NF that runs a nested PF at particular (runtime argument) param values
## paramNodes must be vector of *scalar* nodes
pfLLnf <- nimbleFunction(
    setup = function(model, latentNodes, paramNodes, m) {
        latentNodes <- model$expandNodeNames(latentNodes, sort = TRUE)
        my_PF <- buildPF(model, latentNodes, silent = TRUE)
        if(missing(m)) m <- 10000
    },
    run = function(p = double(1)) {
        values(model, paramNodes) <- p  ## causes ref class '<<-' warning
        ll <- my_PF(m)
        returnType(double())
        return(ll)
    }
)

## specific for psoptim() output: extracts trace information
extractPSOtrace <- function(out) {
    if(is.null(out$stats)) stop("didn't record psoptim() trace")
    x <- t(do.call(cbind, out$stats$x))
    y <- unlist(out$stats$f) * (-1)  ## multiply by -1 to undo fnscale=-1 in psoptim()
    ind <- (y != -Inf)         ## remove -Inf values from y
    y <- y[ind]                ##
    x <- x[ind, , drop=FALSE]  ## (drop=FALSE needed for 1D case!)
    list(x=x, y=y, n=dim(x)[1], p=dim(x)[2])
}

## specific for psoptim() output: 1D or 2D plots of logL estimates
plotPSOtrace <- function(psoTrace, paramNodes) {
    if(psoTrace$p == 1) {
        plot(psoTrace$x[, 1], psoTrace$y, xlab=paramNodes[1], ylab='log-likelihood')
    } else if(psoTrace$p == 2) {
        require(plot3D)
        scatter3D(psoTrace$x[,1], psoTrace$x[,2], psoTrace$y, xlab=paramNodes[1], ylab=paramNodes[2], zlab='log-likelihood')
    }
}




pfLL <- function(...)
    pfLLClass$new(...)



## paramNodes must be vector of *scalar* nodes
pfLLClass <- R6Class(
    'pfLLClass',
    public = list(
        ## class variables
        Rmodel = NULL, RpfLLnf = NULL,
        Cmodel = NULL, CpfLLnf = NULL,
        psoOut = NULL, psoTrace = NULL,
        ## initialize
        initialize = function(
            ## initialize() arguments
            modelarg, latentNodes, paramNodes,
            lower = rep(-Inf,length(paramNodes)), upper = rep(Inf,length(paramNodes)),
            m = 1000, psoIt = 20, setSeed = TRUE, plot = TRUE) {
            ## model and NF creation and compilation
            self$Rmodel <- replicateModel(modelarg)
            self$RpfLLnf <- pfLLnf(self$Rmodel, latentNodes, paramNodes, m)
            self$Cmodel <- compileNimble(self$Rmodel)
            self$CpfLLnf <- compileNimble(self$RpfLLnf, project = self$Rmodel)
            psoControl <- list(fnscale=-1, trace=1, trace.stats=TRUE, REPORT=1, vectorize=TRUE, maxit=psoIt)
            ## pso optimization
            if(setSeed) set.seed(0)
            self$psoOut <- psoptim(paramNodes, self$CpfLLnf$run, lower=lower, upper=upper, control=psoControl)
            self$psoTrace <- extractPSOtrace(self$psoOut)
            if(plot) plotPSOtrace(self$psoTrace, paramNodes)
        }
    )
)





## original version
## good 2D & 3D plotting stuff

## pfLL <- function(model, latent, param, m = 10000, rep = 1, makePlot = TRUE) {
##     pfLLobj <- pfLLClass$new()
##     pfLLobj$initModel(model)
##     pfLLobj$initPF(latent)
##     pfLLobj$initParam(param)
##     pfLLobj$initRep(rep)
##     pfLLobj$calcLL(m)
##     pfLLobj$storeXY()
##     if(makePlot) pfLLobj$makePlot()
##     list(param=pfLLobj$param, ll=pfLLobj$ll, x=pfLLobj$x, y=pfLLobj$y)
## }
## pfLLClass <- R6Class(
##     'pfLLClass',
##     public = list(
##         Rmodel = NULL,
##         Cmodel = NULL,
##         initModel = function(modelarg) {
##             md <- modelarg$modelDef
##             self$Rmodel <- md$newModel()
##             nimCopy(from = modelarg, to = self$Rmodel, logProb = TRUE)
##             for(var in ls(modelarg$isDataEnv)) self$Rmodel$isDataEnv[[var]] <- modelarg$isDataEnv[[var]]
##             self$Cmodel <- compileNimble(self$Rmodel)
##         },
##         Rpf = NULL,
##         Cpf = NULL,
##         initPF = function(latent) {
##             self$Rpf <- buildPF(self$Rmodel, latent)
##             self$Cpf <- compileNimble(self$Rpf, project = self$Rmodel)
##         },
##         param = NULL,
##         paramNames = NULL,
##         nParam = NULL,
##         nParamValues = NULL,
##         initParam = function(param) {
##             self$param <- param
##             self$paramNames <- dimnames(param)[[2]]
##             self$nParam <- length(self$paramNames)
##             self$nParamValues <- dim(param)[1]
##             if(!all(self$paramNames %in% self$Cmodel$getNodeNames(returnScalarComponents = TRUE)))
##                 stop('problem with param names')
##         },
##         rep = NULL,
##         initRep = function(rep) {
##             self$rep <- rep
##         },
##         ll = NULL,
##         calcLL = function(m) {
##             self$ll <- array(NA, c(self$nParamValues, self$rep))
##             for(i in 1:self$nParamValues) {
##                 cat(paste0(round((i-1)/self$nParamValues*100,0), '%\n'))
##                 self$setModelParamValues(row = i)
##                 for(j in 1:self$rep)     self$ll[i, j] <- self$Cpf$run(m)
##             }
##         },
##         setModelParamValues = function(row) {
##             values(self$Cmodel, self$paramNames) <- as.numeric(self$param[row, ])
##         },
##         x = NULL,
##         y = NULL,
##         storeXY = function() {
##             self$x <- array(NA, c(self$nParamValues * self$rep, self$nParam))
##             for(i in 1:self$nParam)   self$x[, i] <- rep(self$param[,i], self$rep)
##             self$y <- as.numeric(self$ll)
##         },
##         makePlot = function() {
##             switch(length(self$paramNames),
##                    plot(self$x[, 1],
##                         self$y,
##                         xlab=self$paramNames,
##                         ylab='log-likelihood estimate'),  ## plot for 1 param
##                    { require(plot3D)
##                      scatter3D(self$x[, 1],
##                                self$x[, 2],
##                                self$y,
##                                xlab=self$paramNames[1],
##                                ylab=self$paramNames[2],
##                                zlab='log-likelihood estimate') }  ## plot for 2 params
##                    )
##         }
##     )
## )











