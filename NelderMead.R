
## seeing what Nelder Mead optimization does with this

myKF <- function(a,b) KF_ll(list(y=data$y, sigOE=0.05, sigPN=0.2, a=a, b=b))
NM <- function(...) NMClass$new(...)
NMClass <- R6Class(
    'NMClass',
    public = list(
        cnf = NULL,
        x = NULL, y = NULL, v = NULL, n = NULL,
        xhist = NULL,
        alpha = NULL, beta = NULL, gamma = NULL, delta = NULL, k = NULL, KF = NULL,
        initialize = function(xinit, cnf, k = 1.25, m = 5000, KF = FALSE) {
            set.seed(0)
            self$cnf <- cnf
            self$k <- k
            self$KF <- KF
            self$cnf$setM(m)
            cnf$reset()
            optim(xinit, cnf$run, control = list(fnscale = -1, maxit = 13))
            ii <- c(1, 5, 6)
            self$x <- cnf$getX()[ii,]
            self$y <- cnf$getY()[ii]
            self$v <- cnf$getV()[ii]
            self$n <- length(self$v)
            self$sortXYV()
            self$xhist <- self$x
            self$alpha <- 1
            self$beta <- 0.5
            self$gamma <- 2
            self$delta <- 0.9##0.5
        },
        sortXYV = function() {
            ix <- sort(self$y, decreasing=TRUE, index.return=TRUE)$ix
            self$x <- self$x[ix,]
            self$y <- self$y[ix]
            self$v <- self$v[ix]
        },
        incPrec = function() {
            message('increasing precision')
            cnf$scaleM(self$k)
            message('m is now: ', cnf$getM())
        },
        decPrec = function() {
            message('decreasing precision')
            cnf$scaleM(1/self$k)
            message('m is now: ', cnf$getM())
        },
        calcOLD = function(xcalc) {
            x <- cnf$getX()
            y <- cnf$getY()
            v <- cnf$getV()
            for(i in 1:dim(x)[1]) {
                if(all(xcalc == x[i,])) {
                    self$xhist <- rbind(self$xhist, xcalc)
                    return(c(y[i], v[i]))
                }
            }
            test <- cnf$run(xcalc)
            if(test == -Inf) return(c(-Inf, Inf))
            print(xcalc)
            print(test)
            stop('coundlt find right xcalc point...')
        },
        calc = function(xcalc) {
            if(self$KF) {
                ll <- myKF(xcalc[1], xcalc[2])
                return(c(ll, 1))
            }
            ll <- cnf$run(xcalc)
            if(ll == -Inf) return(c(-Inf, Inf))
            lastY <- cnf$getLastY(0)
            lastV <- cnf$getLastV(0)
            if(ll != lastY) stop('wrong')
            return(c(lastY, lastV))
        },
        accept = function(xacc, yacc, vacc, pos) {
            self$x[pos,] <- xacc
            self$y[pos] <- yacc
            self$v[pos] <- vacc
            self$xhist <- rbind(self$xhist, xacc)
            self$sortXYV()
        },
        adapt = function() {
            T <- sum(((self$y - mean(self$y)) / sqrt(self$v))^2)
            if(!self$KF) message('T is: ', T) else message('NA')
            pn <- pchisq(T, df=self$n)
            if(!self$KF) message('pn is: ', pn)  else message('NA')
            if(pn < 0.5) {  ## fail to rejust H0: means are equal
                self$incPrec()
                ##if(!self$KF) message('fail to reject H0; would increase precision here') else message('NA')
            } else {
                self$decPrec()
                ##if(!self$KF) message('reject H0; decrease precision here') else message('NA')
            }
        },
        step = function() {
            self$adapt()
            for(i in 1:1) {
                thisx <- self$x[i,]
                thisfv <- self$calc(thisx)
                self$accept(thisx, thisfv[1], thisfv[2], i)
                message('recalc xbest: ', thisfv[1])
                print(thisx)
            }
            self$sortXYV()
            xbest <- self$x[1,]
            xsecond <- self$x[2,]
            y1 <- self$y[1]
            y2 <- self$y[2]
            message('xbest: ', y1)
            print(xbest)
            message('xsecond: ', y2)
            print(xsecond)
            message('xworst: ', self$y[3])
            print(self$x[3,])
            xc <- apply(self$x[-self$n,], 2, mean)
            message('xc:')
            print(xc)
            xworst <- self$x[self$n,]
            yworst <- self$y[self$n]
            xr <- xc + self$alpha*(xc-xworst)
            ##messge('calculate xr')
            fvr <- self$calc(xr)
            fr <- fvr[1]
            message('calc xr: ', fr)
            print(xr)
            if(y1 >= fr & fr > y2) { message('reflection'); self$accept(xr, fr, fvr[2], self$n); return() }
            if(fr > y1) {
                xe <- xc + self$gamma*(xr-xc)
                fve <- self$calc(xe)
                fe <- fve[1]
                message('calc xe: ', fe)
                message('variance of xe calculation: ', fve[2])
                print(xe)
                ##if(fe > fr)  { self$accept(xe, fe, fve[2], self$n); return() }  ## greedy min.
                if(fe > y1)  { message('expansion'); self$accept(xe, fe, fve[2], self$n); return() }  ## greedy exp.
                message('reflection second chance'); self$accept(xr, fr, fvr[2], self$n); return()
            }
            if(fr <= y2) {
                if(fr > yworst) {
                    xcon <- xc + self$beta * (xr - xc)
                    fvcon <- self$calc(xcon)
                    fcon <- fvcon[1]
                    message('calc xcontract 2: ', fcon)
                    print(xcon)
                    ##if(fcon >= fr) { message('contraction 2 from xr'); self$accept(xcon, fcon, fvcon[2], self$n); return() }
                    ## DO NOT ALLOW SHRINKS
                    if(fcon < fr) { message('would have shrunk, not going to') }
                    if(TRUE) { message('contraction 2 from xr'); self$accept(xcon, fcon, fvcon[2], self$n); return() }
                }
                if(fr <= yworst) {
                    xcon <- xc + self$beta * (xworst - xc)
                    fvcon <- self$calc(xcon)
                    fcon <- fvcon[1]
                    message('calc xcontract 1: ', fcon)
                    print(xcon)
                    ##if(fcon > yworst) { message('contraction 1 from xworst'); self$accept(xcon, fcon, fvcon[2], self$n); return() }
                    ## DO NOT ALLOW SHRINKS
                    if(fcon <= yworst) { message('would have shrunk, not going to') }
                    if(TRUE) { message('contraction 1 from xworst'); self$accept(xcon, fcon, fvcon[2], self$n); return() }
                }
            }
            message('shrink')
            for(j in 2:self$n) {
                xnew <- xbest + self$delta * (self$x[j,] - xbest)
                fvnew <- self$calc(xnew)
                fnew <- fvnew[1]
                self$accept(xnew, fnew, fvnew[2], j)
            }
        }
    )
)
