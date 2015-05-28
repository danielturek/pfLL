
require(igraph)
require(nimble)
require(R6)
require(pso)        ## particle swarm optimization psoptim()
##require(mgcv)     ## generalized additive model fitting gam()
require(grDevices)  ## convex hull chull()
require(geometry)   ## multidimensional convex hull and search
library(mvtnorm)    ## multivariate normal generation

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
        d <- length(trans)
        trans <- matrix(as.numeric(c(identity=1, log=2, logit=3)[trans]), c(1,d))
    },
    run = function(vals = double(1), inv = double()) {
        if(length(vals) != d) print('ERROR: wrong length input to trans NF')
        declare(newVals, double(1,d))
        for(i in 1:d) {
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
## param must be vector of *scalar* nodes
pfLLnf <- nimbleFunction(
    setup = function(model, latent, param, trans, d, m) {
        latent <- model$expandNodeNames(latent, sort = TRUE)
        my_trans <- transNF(trans)
        my_PF <- buildPF(model, latent, silent = TRUE)
        if(missing(m)) m <- 10000
        cur <- 0
        max <- increment <- 100  ## must be > 1
        x <- array(0, c(increment, d))
        y <- v <- rep(0, increment)
    },
    run = function(transformedVals = double(1)) {
        vals <- my_trans(transformedVals, 1)
        values(model, param) <<- vals
        ll <- my_PF$run(m)
        if(ll != -Inf) {
            if(cur == max) {
                max <<- max + increment
                xTemp <- x
                yTemp <- y
                vTemp <- v
                setSize(x, max, d)
                setSize(y, max)
                setSize(v, max)
                x[1:cur, 1:d] <<- xTemp
                y[1:cur]      <<- yTemp
                v[1:cur]      <<- vTemp
            }
            cur <<- cur + 1
            x[cur, 1:d] <<- transformedVals
            y[cur]      <<- ll
            v[cur]      <<- my_PF$getVarLL()
        }
        returnType(double())
        return(ll)
    },
    methods = list(
        setM = function(mNew = integer()) m <<- mNew,
        getM = function() { returnType(integer()); return(m)             },
        getX = function() { returnType(double(2)); return(x[1:cur, 1:d]) },
        getY = function() { returnType(double(1)); return(y[1:cur])      },
        getV = function() { returnType(double(1)); return(v[1:cur])      },
        reset = function() {
            cur <<- 0
            max <<- increment
            setSize(x, max, d)
            setSize(y, max)
            setSize(v, max)
            x <<- x * 0
            y <<- y * 0
            v <<- v * 0
        }
    )
)


pfLL <- function(...)
    pfLLClass$new(...)

## param must be vector of *scalar* nodes
pfLLClass <- R6Class(
    'pfLLClass',
    public = list(
        trunc = NULL, plot = NULL,       ## set from initialize()
        param = NULL, d = NULL, numQLMparams = NULL, trans = NULL,            ## set in initParamArgs()
        Rmodel = NULL, rnf = NULL, Cmodel = NULL, cnf = NULL, Ctrans = NULL,  ## set in createNimbleObjects()
        init = NULL, tInit = NULL, tLower = NULL, tUpper = NULL,              ## set in initLowerUpperBest()
        MLE = NULL, tMLE = NULL,
        ## persistent variables:
        bestX = NULL, tBestX = NULL, covX = NULL, bestY = NULL,
        x = NULL, y = NULL, v = NULL, nPoints = NULL,
        quadLM = NULL, quadLMABC = NULL, quadLMcheck = NULL,
        plot2Dtimes = NULL,
        ##
        ## initialization functions
        ##
        initialize = function(
            model, latent, param, lower, upper, init, trans,
            m=5000, trunc=0.9, MLE, setSeed=TRUE, plot=TRUE) {
            if(setSeed) set.seed(0)
            self$trunc <- trunc;     self$plot <- plot
            self$initParamArgs(param, trans)
            self$initNimbleObjects(model, latent, m)
            self$initLowerUpper(init, lower, upper, MLE)
            self$initPersistentVars()
            message('Fitting ', self$d, ' state space model parameters: ', paste0(self$param,collapse=', '))
            message('Requires ', self$numQLMparams, ' parameters for negative-quadratic surface approximation')
        },
        initParamArgs = function(param, trans) {
            self$param <- param;   self$d <- length(self$param);   self$numQLMparams <- self$d*(self$d+3)/2+1
            if(missing(trans)) trans <- rep('identity', self$d)
            if(is.list(trans)) trans <- sapply(trans, function(t) if(identical(t,identity)) 'identity' else if(identical(t,log)) 'log' else if(identical(t,logit)) 'logit' else stop('trans'))
            if(!all(trans %in% c('identity','log','logit'))) stop('trans')
            if(length(trans) != self$d) stop('trans')
            self$trans <- trans
        },
        initNimbleObjects = function(model, latent, m) {
            cat('Compiling model and particle filter...\n')
            self$Rmodel <- replicateModel(model)
            self$rnf <- pfLLnf(self$Rmodel, latent, self$param, self$trans, self$d, m)
            self$Cmodel <- compileNimble(self$Rmodel)
            self$cnf <- compileNimble(self$rnf, project = self$Rmodel)
            self$Ctrans <- self$cnf$my_trans
        },
        initLowerUpper = function(init, lower, upper, MLE) {
            self$init   <- init
            self$tInit  <- self$Ctrans$run(init, 0)
            self$tLower <- self$Ctrans$run(lower, 0)
            self$tUpper <- self$Ctrans$run(upper, 0)
            if(!missing(MLE)) {
                self$MLE <- MLE
                self$tMLE <- self$Ctrans$run(MLE, 0)
            }
        },
        initPersistentVars = function() {
            self$bestX <- self$init;    self$tBestX <- self$Ctrans$run(self$init, 0)
            names(self$bestX) <- self$param;     names(self$tBestX) <- self$param
            self$covX <- diag(self$d)
            self$bestY <- -Inf
            self$resetSamples(print = FALSE)
            self$extractSamples(reduce = FALSE, print = FALSE)    ## sets x, y, v, nPoints
            self$plot2Dtimes <- 0
        },
        ##
        ## internal-only functions
        ##
        internalSetBest = function(newTBestX, newBestY) {
            if(length(newTBestX) != self$d) stop('wrong length for new best')
            newBestX <- self$Ctrans$run(newTBestX, 1)
            names(newTBestX) <- self$param;    names(newBestX) <- self$param
            self$tBestX <- newTBestX;   self$bestX <- newBestX;    self$bestY <- newBestY
        },
        internalSetCovX = function(newCov) {
            if(!all(dim(newCov) == c(self$d, self$d))) stop('wrong dimensions for newCov')
            rownames(newCov) <- colnames(newCov) <- self$param
            self$covX <- newCov
        },
        ##
        ## functions dealing with samples
        ##
        resetSamples = function(print = TRUE) {
            if(print) message('Reset stored samples')
            self$cnf$reset()
        },
        extractSamples = function(reduce = TRUE, print = TRUE) {
            x <- self$cnf$getX();   y <- self$cnf$getY();   v <- self$cnf$getV()
            if(print) message('Extracted ', length(y), ' sample', if(length(y)==1) '' else 's')
            if(reduce) {
                ind <- keepInd <- 2*(max(y)-y) < qchisq(self$trunc, self$d)
                if(sum(ind) > self$d) {
                    xKeep <- x[ind, , drop=FALSE]
                    if(self$d == 1) {
                        keepInd <- (x[,1] >= range(xKeep)[1]) & (x[,1] <= range(xKeep)[2])
                    } else {
                        del <- try(delaunayn(xKeep, options=c('QJ','Pp')), silent = TRUE)
                        if(!inherits(del, 'try-error')) {
                            keepAr <- tsearchn(xKeep, del, x)
                            keepInd <- !is.na(keepAr$idx)
                        }
                    }
                }
                x <- x[keepInd, , drop=FALSE];   y <- y[keepInd];   v <- v[keepInd]
                if(print) message(floor(100*self$trunc), '-percent LRT reduced these to ', length(y), ' sample', if(length(y)==1) '' else 's')
            }
            colnames(x) <- self$param
            self$x <- x;     self$y <- y;     self$v <- v;    self$nPoints <- length(self$y)
        },
        setBestFromSamples = function(print = TRUE) {
            yMax <- max(self$y)
            ind <- which(self$y == yMax)
            txbest <- self$x[ind, ]
            self$internalSetBest(txbest, yMax)
            if(print) { message('Set bestX from samples values:')
                        print(self$bestX)
                        message('Set corresponding best logL: ', self$bestY) }
        },
        setCovXFromSamples = function(print = TRUE) {
            newCov <- cov(self$x)
            self$internalSetCovX(newCov)
            if(print) message('Set covX from sample values')
        },
        ##
        ## functions which evaluate the PF (the only methods which may call self$cnf$run)
        ##
        evaluatePFonXarray = function(xArray, mScale = 1, print = TRUE) {
            origM <- self$cnf$getM();    tempM <- floor(mScale*origM);    self$cnf$setM(tempM)
            if(print) if(mScale != 1) message('Adusting number of particles to: ', tempM)
            numPoints <- dim(xArray)[1];     if(numPoints==0) stop('zero points passed to evalPFonXarray')
            if(print) message('Evaluating particle filter on ', numPoints, ' points...')
            ##for(i in 1:numPoints) self$cnf$run(xArray[i, ])
            lapply(1:numPoints, function(i) self$cnf$run(xArray[i, ]))
            self$cnf$setM(origM);    return(invisible(NULL))
        },
        runOptim = function(numIt, print = TRUE) {
            if(print) message('Running optim for ', numIt, ' iterations...')
            optim(self$tBestX, self$cnf$run, control=list(fnscale=-1, maxit=numIt))
        },
        runPSOptim = function(numIt, print = TRUE) {
            if(print) message('Running PSO for ', numIt, ' iterations...')
            psoptim(self$tBestX, self$cnf$run, lower=self$tLower, upper=self$tUpper,
                    control=list(fnscale=-1, maxit=numIt, vectorize=TRUE))
        },
        ##
        ## fit quadLM
        ##
        fitQuadLM = function(print = TRUE) {
            xterms <- c(self$param, paste0('I(', self$param, '^2)'), if(self$d > 1) combn(self$param, 2, function(x) paste0(x, collapse=':')) else character())
            form <- as.formula(paste0('y ~ ', paste0(xterms, collapse=' + ')))
            df <- as.data.frame(cbind(self$x, y = self$y))
            self$quadLM <- lm(form, weights = 1/self$v, data = df)
            ##if(print) print(summary(self$quadLM))
            if(print) message('Fit ', self$numQLMparams, '-parameter negative-quadratic linear model to ', self$nPoints, ' data points')
            if(print) message('Checking model fit')
            self$quadLMcheck <- 1
            self$quadLMABC <- NULL
            coef <- self$quadLM$coef
            A <- array(0, c(self$d, self$d), dimnames = list(self$param, self$param))
            for(i in seq_along(self$param)) {
                A[i, i] <- coef[paste0('I(', self$param[i], '^2)')]
                if(i!=self$d) for(j in (i+1):self$d) A[i, j] <- A[j, i] <- 1/2 * coef[paste0(self$param[i], ':', self$param[j])]
            }
            b <- array(coef[self$param], c(self$d, 1), dimnames = list(self$param, NULL))
            c <- coef['(Intercept)']
            if(any(is.na(coef))) { if(print) message('Model coefficient is NA: singular design matrix'); self$quadLMcheck <- -1 }
            if(any(A==0) || any(b==0) || (c==0) ) stop('Coefficients not extracted correctly')
            if(any(eigen(A)$values > 0)) { if(print) message('Negative-quadratic model fit is NOT strictly concave down'); self$quadLMcheck <- -2 }
            self$quadLMABC <- list(A=A, b=b, c=c)
            if(print) if(self$quadLMcheck == 1) message('Model fit OK')
        },
        setBestXCovXFromQuadLM = function(print = TRUE) {
            A <- self$quadLMABC$A;   b <- self$quadLMABC$b;   c <- self$quadLMABC$c
            paramT <- t(-1/2 * solve(A) %*% b)[1, ]
            logL <- (-1/4 * t(b) %*% solve(A) %*% b + c)[1,1]
            self$internalSetBest(paramT, logL)
            newCov <- -1/2 * solve(A)
            self$internalSetCovX(newCov)
            if(print) { message('Set bestX from negative-quadratic linear model')
                        print(self$bestX)
                        message('Set corresponding best logL: ', self$bestY)
                        message('Set covX from quadratic linear model fit') }
        },
        ##
        ## loop optim
        ##
        loopOptim = function(iter = 10, mult = 1, alg = 'optim', print = TRUE) {
            minPoints <- floor(self$numQLMparams * mult)
            if(!(alg %in% c('optim', 'PSO'))) stop('unknown algorithm')
            if(print) message('Loop ', alg, ' until we collect ', minPoints, ' points within ', floor(100*self$trunc), '-percent LRT interval')
            thisIter <- 0
            while(self$nPoints < minPoints) {
                thisIter <- thisIter + iter
                if(alg=='optim')   self$runOptim(thisIter, print = print)
                if(alg=='PSO')     self$runPSOptim(thisIter, print = print)
                self$extractSamples(print = print)
                self$setBestFromSamples(print = FALSE)
                if(print) { message('Best parameter location found thus far:')
                            print(self$bestX)
                            message('Best logL found there: ', self$bestY) }
            }
            self$setCovXFromSamples(print = print)
            if(self$plot) self$plot2D()
        },
        ##
        ## MVN cloud
        ##
        mvnCloud = function(pointsScale = 10, mScale = 1, print = TRUE) {
            numPoints <- self$numQLMparams * pointsScale
            if(print) message('Generating MVN point cloud of ', numPoints, ' points')
            xNew <- rmvnorm(numPoints, self$tBestX, self$covX)
            if(print) message('Evaluating particle filter on point cloud...')
            self$evaluatePFonXarray(xNew, mScale = mScale, print = FALSE)
            self$extractSamples(print = print)
            self$setBestFromSamples(print = FALSE)
            if(self$plot) self$plot2D(newPlot = FALSE, bestDot = FALSE)
        },
        ##
        ## linear grid in directions of positive eigen values
        ##
        linearGrid = function(numPoints = 100, print = TRUE) {
            if(print) message('Expanding linear search grid in direction of each positive eigenvalue')
            eigA <- eigen(self$quadLMABC$A)
            evMat <- eigA$vectors
            reParamCovX <- evMat %*% self$covX %*% t(evMat)   ## this way is correct
            posInd <- which(eigA$values > 0);   if(length(posInd)==0) stop('')
            if(print) message('Identified ', length(posInd), ' increasing direction', if(length(posInd)>1) 's' else '')
            numPoints <- 100
            for(iPosInd in seq_along(posInd)) {
                ei <- posInd[iPosInd]
                linPoints <- matrix(1:numPoints) * sqrt(reParamCovX[ei,ei])/10
                xNew <- linPoints %*% evMat[,ei] + matrix(1,numPoints) %*% self$tBestX
                if(self$plot) if(dev.cur()!=1) points(xNew[,1], xNew[,2], col='purple', pch=18)
                if(print) message('Direction ', iPosInd, ': evaluating particle filter on ', numPoints, ' new points...')
                self$evaluatePFonXarray(xNew, mScale = 1, print = FALSE)
            }
        },
        plot2D = function(x, newPlot=TRUE, initCol='red', mleCol='green', bestCol='blue', bestDot=TRUE) {
            self$plot2Dtimes <- self$plot2Dtimes + 1
            colInd <- (self$plot2Dtimes %% 4) + 5  ## cycles through colors 5,6,7,8
            if(missing(x))  x <- self$x
            if(newPlot) {
                rng <- apply(x, 2, range)
                xlim <- c(rng[1,1], rng[2,1])
                ylim <- c(rng[1,2], rng[2,2])
                if(!is.null(self$MLE)) {
                    warning('change the bounds so it over-extends past the MLE')
                    xlim <- c(min(xlim[1], self$tMLE[1]), max(xlim[2], self$tMLE[1]))
                    ylim <- c(min(ylim[1], self$tMLE[2]), max(ylim[2], self$tMLE[2]))
                }
                dev.new()
                plot(x[,1], x[,2], col=colInd, xlab=self$param[1], ylab=self$param[2], xlim=xlim, ylim=ylim)
            }
            points(x[,1], x[,2], col=colInd)
            self$addDotAtInit(initCol)
            self$addDotAtMLE(mleCol)
            if(bestDot) self$addDotAtBest(bestCol)
        },
        addDotAtInit = function(initCol = 'red') {
            points(self$tInit[1],  self$tInit[2],  pch=19, col=initCol)
        },
        addDotAtMLE = function(mleCol = 'green') {
            if(!is.null(MLE))
                points(self$tMLE[1], self$tMLE[2], pch=19, col=mleCol)
        },
        addDotAtBest = function(bestCol = 'blue') {
            points(self$tBestX[1], self$tBestX[2], pch=19, col=bestCol)
        },
        plotProjections = function() {
            message('Generating projection plots...')
            yPred <- as.numeric(predict(self$quadLM))
            numCol <- if(self$d == 1) 1 else if(self$d < 7) 2 else if(self$d < 13) 3 else 4
            numRow <- floor((self$d-1)/numCol) + 1
            dev.new(); par(mfrow = c(numRow, numCol))
            for(iParam in 1:self$d) {
                xlab <- if(self$trans[iParam]=='identity') self$param[iParam] else paste0(self$trans[iParam],'(',self$param[iParam],')')
                ylim <- c(min(self$y), max(c(self$y, yPred)))
                plot(self$x[,iParam], self$y, xlab=xlab, ylab='log-likelihood', ylim=ylim, pch=20)
                segments(x0=self$x[,iParam], y0=self$y-2*sqrt(self$v), y1=self$y+2*sqrt(self$v), col='grey64')
                points(self$x[,iParam], self$y, pch=20)
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
            self$Rmodel <- self$Cmodel <- self$rnf <- self$cnf <- self$Ctrans <- NULL
        }
        ## if(model == 'GAM') {
        ##     warning('GAM model fitting is in total dis-repair')
        ##     kMax <- min(floor(length(self$y)^(1/self$d)), 6)
        ##     if(kMax < 3) stop('Need more data to fit GAM; increase psoIt argument')
        ##     form <- as.formula(paste0('y ~ te(', paste0(self$param, collapse=', '), ', k = ', kMax,')'))
        ##     df <- as.data.frame(cbind(self$x, y = self$y))
        ##     self$fittedModel <- gam(form, weights = 1/self$v, data = df)
        ##     if(self$d == 2)     { dev.new();   plot(self$fittedModel) }
        ## }
        ## optimFxn = function(vals) {
        ##     if(self$d == 1) { if((vals < self$cHull[1]) | (vals > self$cHull[2])) return(-Inf)
        ##                   } else { if(is.na(tsearchn(self$cHull[[1]], self$cHull[[2]], t(matrix(vals)))$idx)) return(-Inf) }
        ##     df <- data.frame(as.list(vals))
        ##     names(df) <- self$param
        ##     as.numeric(predict(self$fittedModel, df))
        ## },
        ## runSurfaceOptim = function() {
        ##     cat('Optimizing to find fitted-surface maximum...\n')
        ##     if(is.null(self$cHull)) stop('never assigned self$cHull; this should never happen')
        ##     if(self$d == 1) {
        ##         self$optimOut <- optimize(self$optimFxn, range(self$x), maximum=TRUE)
        ##         paramT <- self$optimOut$maximum; logL <- self$optimOut$objective
        ##     } else {
        ##         self$optimOut <- optim(self$bestX, self$optimFxn, control=list(fnscale=-1))
        ##         if(self$optimOut$convergence != 0)   warning('optim() failed to converge')
        ##         paramT <- self$optimOut$par;     logL <- self$optimOut$value
        ##     }
        ##     param <- self$Ctrans$run(paramT, 1)
        ##     names(param) <- self$param;   names(paramT) <- paste0(self$param, 'T')
        ##     self$max <- list(param=param, paramT=paramT, logL=logL)
        ##     cat('fitted surface peak location:\n');       print(self$max$param)
        ##     cat(paste0('fitted surface peak value: ', self$max$logL, '\n'))
        ## },
    )
)


## NOTES:
## removing -Inf values from y
## truncating y 99% profile-likelihood CI
## then also keeping x & y pairs inside the resulting p-dim convex hull
## estimating SE of PF predictions, using as 'weights' for lm()
## restricting model predictions to inside the same p-dim convex hull

## deprecated:
## truncating y to top 80%
## kMax for gam() set to 4 to avoid over-fitting
## restrict gam() model predictions to within range of x

KF_ll <- function(lst) {
    y <- lst$y
    t <- length(y)
    b <- lst$b
    if(!is.null(lst$a)) {
        a <- lst$a
    } else if(!is.null(lst$mu)) {
        a <- 1 - b/lst$mu
    } else stop('lst missing both a and mu')
    sigPN <- lst$sigPN
    sigPN2 <- sigPN*sigPN
    sigOE <- lst$sigOE
    sigOE2 <- sigOE*sigOE
    mu_x <- b/(1-a)
    var_x <- sigPN2 / (1-a^2)
    cov_xy <- var_x
    var_y <- var_x + sigOE2
    ll <- dnorm(mu_x, y[1], sqrt(var_y), log = TRUE)
    ##message('iNode: 1,     logL: ', ll)
    mu_x <- mu_x + (cov_xy / var_y) * (y[1] - mu_x)
    var_x <- var_x - cov_xy*cov_xy / var_y
    for(i in 2:t) {
        mu_x <- a*mu_x + b
        var_x <- a^2 * var_x + sigPN2
        if(!is.na(y[i])) {
            cov_xy <- var_x
            var_y <- var_x + sigOE2
            ll <- ll + dnorm(mu_x, y[i], sqrt(var_y), log = TRUE)
            ##if(i <= 3)  message('iNode: ', i, ',     logL: ', dnorm(mu_x, y[i], sqrt(var_y), log = TRUE), ',     cumulative: ', ll)
            mu_x <- mu_x + (cov_xy / var_y) * (y[i] - mu_x)
            var_x <- var_x - cov_xy*cov_xy / var_y
        }
    }
    ll
}








