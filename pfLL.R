
require(igraph)
require(nimble)
require(R6)
require(pso)        ## particle swarm optimization psoptim()
require(mgcv)       ## generalized additive model fitting gam()
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

## parameter transformations and inverse transformations
## 1:identity, 2:log, 3:logit
transNF <- nimbleFunction(
    setup = function(trans) {
        p <- length(trans)
        trans <- matrix(as.numeric(c(identity=1, log=2, logit=3)[trans]), c(1,p))
    },
    run = function(vals = double(1), inv = double()) {
        if(length(vals) != p) print('ERROR: wrong length input to trans NF')
        declare(newVals, double(1,p))
        for(i in 1:p) {
            if(trans[1,i] == 1) {        ## identity
                newVals[i] <- vals[i]
            } else if(trans[1,i] == 2) { ## log
                if(inv==0) newVals[i] <- log(vals[i])
                if(inv==1) newVals[i] <- exp(vals[i])
            } else if(trans[1,i] == 3) { ## logit
                if(inv==0) newVals[i] <- logit(vals[i])
                if(inv==1) newVals[i] <- expit(vals[i])
            }
        }
        returnType(double(1))
        return(newVals)
    }
)

## runs a nested PF at particular (runtime argument) param values
## runtime argument of param values are on *transformed* scales
## paramNodes must be vector of *scalar* nodes
pfLLnf <- nimbleFunction(
    setup = function(model, latentNodes, paramNodes, trans, m) {
        latentNodes <- model$expandNodeNames(latentNodes, sort = TRUE)
        my_trans <- transNF(trans)
        my_PF <- buildPF(model, latentNodes, silent = TRUE)
        if(missing(m)) m <- 10000
    },
    run = function(transformedVals = double(1)) {
        vals <- my_trans(transformedVals, 1)
        values(model, paramNodes) <<- vals
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
        paramNodes = NULL, p = NULL, trans = NULL,
        Rmodel = NULL, RpfLLnf = NULL, Cmodel = NULL, CpfLLnf = NULL, Ctrans = NULL,
        psoOut = NULL, x = NULL, y = NULL,
        ## quadLM = NULL,
        GAM = NULL,
        optimOut = NULL, max = NULL,
        initialize = function(
            modelarg, latentNodes, paramNodes,
            lower = rep(-Inf, length(paramNodes)), upper = rep(Inf, length(paramNodes)), trans,
            m = 1000, psoIt = 20, trunc = 0.8, setSeed = TRUE, plot = TRUE) {
            self$processParamArgs(paramNodes, trans)
            self$createNimbleObjects(modelarg, latentNodes, m)
            self$runPSOptim(lower, upper, psoIt, setSeed)
            self$extractTrace(trunc)
            ##self$fitQuadraticLM()
            self$fitGAM()
            self$runSurfaceOptim()
            if(plot) self$plot()
            self$cleanup()
        },
        processParamArgs = function(paramNodes, trans) {
            self$paramNodes <- paramNodes
            self$p <- length(self$paramNodes)
            if(missing(trans)) trans <- rep('identity', self$p)
            if(is.list(trans)) trans <- sapply(trans, function(t) if(identical(t,identity)) 'identity' else if(identical(t,log)) 'log' else if(identical(t,logit)) 'logit' else stop('trans'))
            if(!all(trans %in% c('identity','log','logit'))) stop('trans')
            if(length(trans) != self$p) stop('trans')
            self$trans <- trans
        },
        createNimbleObjects = function(modelarg, latentNodes, m) {
            self$Rmodel <- replicateModel(modelarg)
            self$RpfLLnf <- pfLLnf(self$Rmodel, latentNodes, self$paramNodes, self$trans, m)
            self$Cmodel <- compileNimble(self$Rmodel)
            self$CpfLLnf <- compileNimble(self$RpfLLnf, project = self$Rmodel)
            self$Ctrans <- self$CpfLLnf$my_trans
        },
        runPSOptim = function(lower, upper, psoIt, setSeed) {
            tLower <- self$Ctrans$run(lower, 0);   tUpper <- self$Ctrans$run(upper, 0)
            psoControl <- list(fnscale=-1, trace=1, trace.stats=TRUE, REPORT=1, vectorize=TRUE, maxit=psoIt)
            if(setSeed) set.seed(0)
            self$psoOut <- psoptim(self$paramNodes, self$CpfLLnf$run, lower=tLower, upper=tUpper, control=psoControl)
        },
        extractTrace = function(trunc) {
            x <- t(do.call(cbind, self$psoOut$stats$x))
            colnames(x) <- self$paramNodes
            y <- unlist(self$psoOut$stats$f) * (-1)  ## (-1) to undo fnscale=-1 in psoptim()
            x <- x[(y != -Inf), , drop=FALSE];    y <- y[(y != -Inf)]   ## remove -Inf from y
            if(trunc < 0 | trunc > 1) stop('bad value of trunc')
            minY <- quantile(y, 1-trunc)
            x <- x[(y > minY), , drop = FALSE];   y <- y[(y > minY)]
            self$x <- x;                          self$y <- y
        },
        ## fitQuadraticLM = function() {
        ##     xterms <- c(self$paramNodes, paste0('I(', self$paramNodes, '^2)'), if(self$p > 1) combn(self$paramNodes, 2, function(x) paste0(x, collapse=':')) else character())
        ##     form <- as.formula(paste0('y ~ ', paste0(xterms, collapse=' + ')))
        ##     df <- as.data.frame(cbind(self$x, y = self$y))
        ##     self$quadLM <- lm(form, data = df)
        ##     print(summary(self$quadLM))
        ## },
        fitGAM = function() {
            kMax <- min(floor(length(self$y)^(1/self$p)), 6)
            if(kMax < 3) stop('Need more data to fit GAM; increase psoIt argument')
            form <- as.formula(paste0('y ~ te(', paste0(self$paramNodes, collapse=', '), ', k = ', kMax,')'))
            self$GAM <- gam(form, data = as.data.frame(cbind(self$x, y=self$y)))
            if(self$p == 2)     { dev.new();   plot(self$GAM) }
        },
        optimFxn = function(vals) {
            ## IMPORTANT: restrict fitted model predictions to within the range of x
            rng <- apply(self$x, 2, range)
            if(any(vals < rng[1, ])) return(-Inf)
            if(any(vals > rng[2, ])) return(-Inf)
            df <- data.frame(as.list(vals))
            names(df) <- self$paramNodes
            ##as.numeric(predict(self$quadLM, df))
            as.numeric(predict(self$GAM, df))
        },
        runSurfaceOptim = function() {
            if(self$p == 1) {
                self$optimOut <- optimize(self$optimFxn, range(self$x), maximum=TRUE)
                paramT <- self$optimOut$maximum; logL <- self$optimOut$objective
            } else {
                self$optimOut <- optim(self$psoOut$par, self$optimFxn, control=list(fnscale=-1))
                if(self$optimOut$convergence != 0)   warning('optim() failed to converge')
                paramT <- self$optimOut$par;     logL <- self$optimOut$value
            }
            param <- self$Ctrans$run(paramT, 1)
            names(param) <- self$paramNodes;   names(paramT) <- paste0(self$paramNodes, 'T')
            self$max <- list(param=param, paramT=paramT, logL=logL)
            cat('Fitted surface peak location:\n');       print(self$max$param)
            cat('Fitted surface peak value (log-likelihood):\n');   print(self$max$logL)
        },
        plot = function() {
            ##yPred <- as.numeric(predict(self$quadLM))
            yPred <- as.numeric(predict(self$GAM))
            numCol <- if(self$p == 1) 1 else if(self$p < 7) 2 else if(self$p < 13) 3 else 4
            numRow <- floor((self$p-1)/numCol) + 1
            dev.new(); par(mfrow = c(numRow, numCol))
            for(iParam in 1:self$p) {
                xlab <- if(self$trans[iParam]=='identity') self$paramNodes[iParam] else paste0(self$trans[iParam],'(',self$paramNodes[iParam],')')
                ylim <- c(min(self$y), max(c(self$y, yPred)))
                plot(self$x[,iParam], self$y, xlab=xlab, ylab='log-likelihood', ylim=ylim, pch=20)
                sortedX <- sort(self$x[,iParam], index.return = TRUE)
                thisX <- sortedX$x; thisY <- yPred[sortedX$ix]
                points(thisX, thisY, col='blue', pch=4)
                chInd <- chull(thisX, thisY)
                while(chInd[1]!=min(chInd)) chInd <- chInd[c(2:length(chInd), 1)]
                chInd <- chInd[1:which(chInd==max(chInd))]
                lines(thisX[chInd], thisY[chInd], col='red', lwd=2)
                points(x=self$max$paramT[iParam], y=self$max$logL, col='green', pch=19)
            }
        },
        cleanup = function() {
            self$Rmodel <- self$Cmodel <- self$RpfLLnf <- self$CpfLLnf <- self$Ctrans <- NULL
        }
    )
)


## NOTES:
## removing -Inf values from y
## truncating y to 80%
## kMax for gam() set to 4 to avoid over-fitting
## restrict gam() model predictions to within range of x

KF_ll <- function(d) {
    mu <- d$mu
    a <- d$a
    b <- mu*(1-a)
    sigPN <- d$sigPN
    sigPN2 <- sigPN*sigPN
    sigOE <- d$sigOE
    sigOE2 <- sigOE*sigOE
    y <- d$y
    t <- d$t
    mu_x <- b/(1-a)
    var_x <- sigPN2 / (1-a^2)
    cov_xy <- var_x
    var_y <- var_x + sigOE2
    ll <- dnorm(mu_x, y[1], sqrt(var_y), log = TRUE)
    mu_x <- mu_x + (cov_xy / var_y) * (y[1] - mu_x)
    var_x <- var_x - cov_xy*cov_xy / var_y
    for(i in 2:t) {
        mu_x <- a*mu_x + b
        var_x <- a^2 * var_x + sigPN2
        if(!is.na(y[i])) {
            cov_xy <- var_x
            var_y <- var_x + sigOE2
            ll <- ll + dnorm(mu_x, y[i], sqrt(var_y), log = TRUE)
            mu_x <- mu_x + (cov_xy / var_y) * (y[i] - mu_x)
            var_x <- var_x - cov_xy*cov_xy / var_y
        }
    }
    ll
}



