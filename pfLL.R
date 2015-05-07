
require(nimble)
require(R6)
require(pso)  ## particle swarm optimization psoptim()
library(grDevices)  ## convex hull chull()

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
        values(model, paramNodes) <<- p
        ll <- my_PF(m)
        returnType(double())
        return(ll)
    }
)



pfLL <- function(...)
    pfLLClass$new(...)

## paramNodes must be vector of *scalar* nodes
pfLLClass <- R6Class(
    'pfLLClass',
    public = list(
        paramNodes = NULL, p = NULL,
        Rmodel = NULL, RpfLLnf = NULL,
        Cmodel = NULL, CpfLLnf = NULL,
        psoOut = NULL, x = NULL, y = NULL,
        quadLM = NULL,
        initialize = function(
            modelarg, latentNodes, paramNodes,
            lower = rep(-Inf,length(paramNodes)), upper = rep(Inf,length(paramNodes)),
            m = 1000, psoIt = 20, trunc = 0.8, setSeed = TRUE, plot = TRUE) {
            self$paramNodes <- paramNodes
            self$p <- length(self$paramNodes)
            self$Rmodel <- replicateModel(modelarg)
            self$RpfLLnf <- pfLLnf(self$Rmodel, latentNodes, self$paramNodes, m)
            self$Cmodel <- compileNimble(self$Rmodel)
            self$CpfLLnf <- compileNimble(self$RpfLLnf, project = self$Rmodel)
            psoControl <- list(fnscale=-1, trace=1, trace.stats=TRUE, REPORT=1, vectorize=TRUE, maxit=psoIt)
            if(setSeed) set.seed(0)
            self$psoOut <- psoptim(self$paramNodes, self$CpfLLnf$run, lower=lower, upper=upper, control=psoControl)
            self$extractPSOtrace(trunc)
            self$fitQuadraticLM()
            print(summary(self$quadLM))
            if(plot) self$plotPSOtrace()  ## plot 1D, 2D cases only
            self$Rmodel <- self$Cmodel <- self$RpfLLnf <- self$CpfLLnf <- NULL  ## cleanup
        },
        extractPSOtrace = function(trunc) {
            if(is.null(self$psoOut$stats)) stop("didn't record psoptim() trace")
            x <- t(do.call(cbind, self$psoOut$stats$x))
            colnames(x) <- self$paramNodes
            y <- unlist(self$psoOut$stats$f) * (-1)  ## (-1) to undo fnscale=-1 in psoptim()
            x <- x[(y != -Inf), , drop=FALSE]  ## drop=FALSE needed for 1D case
            y <- y[(y != -Inf)]                ## remove -Inf values from y
            if(trunc < 0 | trunc > 1) stop('bad value of trunc')
            minY <- quantile(y, 1-trunc)
            x <- x[(y > minY), , drop = FALSE]
            y <- y[(y > minY)]
            self$x <- x;     self$y <- y
        },
        fitQuadraticLM = function() {
            xs <- paste0('x[, ', 1:self$p, ']')
            xterms <- c(self$paramNodes,
                        paste0('I(', self$paramNodes, '^2)'),
                        if(self$p > 1) combn(self$paramNodes, 2, function(x) paste0(x, collapse=':')) else character()
                        )
            form <- as.formula(paste0('y ~ ', paste0(xterms, collapse=' + ')))
            df <- as.data.frame(cbind(self$x, y = self$y))
            self$quadLM <- lm(form, data = df)
        },
        plotPSOtrace = function() {
            yPred <- as.numeric(predict(self$quadLM))
            numCol <- if(self$p == 1) 1 else if(self$p < 7) 2 else if(self$p < 13) 3 else 4
            numRow <- floor((self$p-1)/numCol) + 1
            par(mfrow = c(numRow, numCol))
            for(iParam in 1:self$p) {
                plot(self$x[,iParam], self$y, xlab=self$paramNodes[iParam], ylab='log-likelihood')
                sortedX <- sort(self$x[,iParam], index.return = TRUE)
                thisX <- sortedX$x; thisY <- yPred[sortedX$ix]
                points(thisX, thisY, col='blue', pch=4)
                hullInd <- chull(thisX, thisY)
                while(hullInd[1]!=min(hullInd)) hullInd <- c(hullInd[2:length(hullInd)], hullInd[1])
                hullInd <- hullInd[1:which(hullInd==max(hullInd))]
                hullX <- thisX[hullInd]; hullY <- thisY[hullInd]
                lines(hullX, hullY, col='red', lwd=2)
            }
        }
    )
)




## second version of the plotting code
##
## if(self$p == 1) {
##     plot(self$x[, 1], self$y, xlab=self$paramNodes[1], ylab='log-likelihood')
##     sortedX <- sort(self$x[, 1], index.return = TRUE)
##     lines(sortedX$x, predict(self$quadLM)[sortedX$ix], col='red')
## } else if(self$p == 2) {
##     require(plot3D)
##     scatter3D(self$psoTrace$x[,1], self$psoTrace$x[,2], self$psoTrace$y, xlab=self$paramNodes[1], ylab=self$paramNodes[2], zlab='log-likelihood')
##     scatter3D(self$psoTrace$x[,1], self$psoTrace$x[,2], predict(self$quadLM), xlab=self$paramNodes[1], ylab=self$paramNodes[2], zlab='log-likelihood', add=TRUE, pch=18)
## }



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











