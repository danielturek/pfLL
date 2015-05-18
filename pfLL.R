
require(igraph)
require(nimble)
require(R6)
require(pso)        ## particle swarm optimization psoptim()
##require(mgcv)       ## generalized additive model fitting gam()
require(grDevices)  ## convex hull chull()
require(geometry)   ## multidimensional convex hull and search

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
        x <- array(as.numeric(NA), c(increment, d))
        y <- v <- rep(as.numeric(NA), increment)
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
        getX = function() { returnType(double(2)); return(x[1:cur, 1:d]) },
        getY = function() { returnType(double(1)); return(y[1:cur])      },
        getV = function() { returnType(double(1)); return(v[1:cur])      }
    )
)


pfLL <- function(...)
    pfLLClass$new(...)

## param must be vector of *scalar* nodes
pfLLClass <- R6Class(
    'pfLLClass',
    public = list(
        param = NULL, d = NULL, trans = NULL,
        Rmodel = NULL, RpfLLnf = NULL, Cmodel = NULL, CpfLLnf = NULL, Ctrans = NULL,
        tLower = NULL, tUpper = NULL, psoIt = NULL, trunc = NULL,
        x = NULL, y = NULL, v = NULL, bestX = NULL, nPoints = NULL,
        A <- NULL, fittedModel = NULL,  max = NULL,
        initialize = function(
            modelarg, latent, param,
            lower = rep(-Inf, length(param)), upper = rep(Inf, length(param)), trans,
            m = 5000, psoIt = 10, trunc = 0.999, setSeed = TRUE, plot = TRUE) {
            if(setSeed) set.seed(0)
            self$initParamArgs(param, trans)
            self$createNimbleObjects(modelarg, latent, m)
            self$initPSArgs(lower, upper, psoIt, trunc)
            while(self$nPoints < 100*self$d) { self$runPSOptim(); self$extractXYV() }
            self$fitModel('quadLM')
            ##self$runSurfaceOptim()
            if(plot) self$plot()
            ##self$cleanup()
        },
        initParamArgs = function(param, trans) {
            self$param <- param;   self$d <- length(self$param)
            if(missing(trans)) trans <- rep('identity', self$d)
            if(is.list(trans)) trans <- sapply(trans, function(t) if(identical(t,identity)) 'identity' else if(identical(t,log)) 'log' else if(identical(t,logit)) 'logit' else stop('trans'))
            if(!all(trans %in% c('identity','log','logit'))) stop('trans')
            if(length(trans) != self$d) stop('trans')
            self$trans <- trans
        },
        createNimbleObjects = function(modelarg, latent, m) {
            cat('Compiling model and particle filter algorithm...\n')
            self$Rmodel <- replicateModel(modelarg)
            self$RpfLLnf <- pfLLnf(self$Rmodel, latent, self$param, self$trans, self$d, m)
            self$Cmodel <- compileNimble(self$Rmodel)
            self$CpfLLnf <- compileNimble(self$RpfLLnf, project = self$Rmodel)
            self$Ctrans <- self$CpfLLnf$my_trans
        },
        initPSArgs = function(lower, upper, psoIt, trunc) {
            self$tLower <- self$Ctrans$run(lower, 0);   self$tUpper <- self$Ctrans$run(upper, 0)
            self$psoIt  <- psoIt;   self$trunc <- trunc
            self$bestX <- rep(NA, self$d);   self$nPoints <- 0
        },
        runPSOptim = function() {
            cat(paste0('Running particle swarm optimizer for ', self$psoIt, ' iterations...\n'))
            psoptim(self$bestX, self$CpfLLnf$run, lower=self$tLower, upper=self$tUpper,
                    control=list(fnscale=-1, vectorize=TRUE, maxit=self$psoIt))
            self$psoIt <- self$psoIt * 2
        },
        extractXYV = function() {
            x <- self$CpfLLnf$getX();   y <- self$CpfLLnf$getY();   v <- self$CpfLLnf$getV()
            cat(paste0('Collected a total of ', length(y), ' finite points, '))
            colnames(x) <- self$param
            ind <- 2*(max(y)-y) < qchisq(self$trunc, self$d)  ## indicies of logL within trunc CI
            if(sum(ind) <= self$d) {
                cat('too few to create convex hull, ');     keepInd <- ind
            } else {
                xKeep <- x[ind,,drop=FALSE]    ## x values corresponding to these logL's
                if(self$d == 1) {
                    keepInd <- (x[,1] >= range(xKeep)[1]) & (x[,1] <= range(xKeep)[2])
                } else {
                    del <- try(delaunayn(xKeep), silent = TRUE)
                    if(inherits(del, 'try-error')) {
                        message('co-planar points, unable to create convex hull, ');     keepInd <- ind
                    } else {
                        keepAr <- tsearchn(xKeep, del, x)
                        keepInd <- !is.na(keepAr$idx)
                    }
                }
            }
            x <- x[keepInd,,drop=FALSE];   y <- y[keepInd];   v <- v[keepInd]
            ## important: not doing next line, too restictive, and fails in co-planar case
            ##self$tLower <- apply(x, 2, min);     self$tUpper <- apply(x, 2, max)
            self$bestX <- x[which(y==max(y)), ];   self$nPoints <- length(y)
            cat(paste0('keeping ', self$nPoints, '.\n'))
            self$x <- x;     self$y <- y;     self$v <- v
        },
        fitModel = function(model) {
            if(model == 'quadLM') {
                cat(paste0('Fitting negative-quadratic model to ', self$nPoints, ' data points...\n'))
                xterms <- c(self$param, paste0('I(', self$param, '^2)'), if(self$d > 1) combn(self$param, 2, function(x) paste0(x, collapse=':')) else character())
                form <- as.formula(paste0('y ~ ', paste0(xterms, collapse=' + ')))
                df <- as.data.frame(cbind(self$x, y = self$y))
                self$fittedModel <- lm(form, weights = 1/self$v, data = df)
                coef <- self$fittedModel$coef
                A <- array(0, c(self$d, self$d), dimnames = list(self$param, self$param))
                for(i in seq_along(self$param)) {
                    A[i, i] <- coef[paste0('I(', self$param[i], '^2)')]
                    if(i!=self$d) for(j in (i+1):self$d) A[i, j] <- A[j, i] <- 1/2 * coef[paste0(self$param[i], ':', self$param[j])]
                }
                b <- array(coef[self$param], c(self$d, 1), dimnames = list(self$param, NULL))
                c <- coef['(Intercept)']
                if(any(eigen(A)$values > 0)) warning('Negative-quadratic model fit is NOT strictly concave down')
                paramT <- t(-1/2 * solve(A) %*% b)[1, ]
                param <- self$Ctrans$run(paramT, 1);    names(param) <- self$param
                logL <- (-1/4 * t(b) %*% solve(A) %*% b + c)[1,1]
                self$A <- A
                self$max <- list(param=param, paramT=paramT, logL=logL)
                print(summary(self$fittedModel))
                cat('fitted surface peak location:\n');   print(self$max$param)
                cat(paste0('fitted surface peak value: ', self$max$logL, '\n'))
            } else if(model == 'GAM') {
                warning('GAM model fitting is in total dis-repair')
                kMax <- min(floor(length(self$y)^(1/self$d)), 6)
                if(kMax < 3) stop('Need more data to fit GAM; increase psoIt argument')
                form <- as.formula(paste0('y ~ te(', paste0(self$param, collapse=', '), ', k = ', kMax,')'))
                df <- as.data.frame(cbind(self$x, y = self$y))
                self$fittedModel <- gam(form, weights = 1/self$v, data = df)
                if(self$d == 2)     { dev.new();   plot(self$fittedModel) }
            } else stop('unknown model')
        },
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
        plot = function() {
            cat('Generating projection plots...\n')
            yPred <- as.numeric(predict(self$fittedModel))
            numCol <- if(self$d == 1) 1 else if(self$d < 7) 2 else if(self$d < 13) 3 else 4
            numRow <- floor((self$d-1)/numCol) + 1
            dev.new(); par(mfrow = c(numRow, numCol))
            for(iParam in 1:self$d) {
                xlab <- if(self$trans[iParam]=='identity') self$param[iParam] else paste0(self$trans[iParam],'(',self$param[iParam],')')
                ylim <- c(min(self$y), max(c(self$y, yPred)))
                plot(self$x[,iParam], self$y, xlab=xlab, ylab='log-likelihood', ylim=ylim, pch=20)
                segments(x0=self$x[,iParam], y0=self$y-sqrt(self$v)/2, y1=self$y+sqrt(self$v)/2)
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



