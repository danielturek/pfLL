
## basis for natural cubic spline on [0, 1]
basis <- function(x, z) {
    ((x-1/2)^2 - 1/12) * ((z-1/2)^2 - 1/12) / 4 -
        ((abs(x-z) - 1/2)^4 - 1/2*(abs(x-z) - 1/2)^2 + 7/240) / 24
}

splineX <- function(xs, xknot) {
    q <- length(xknot) + 2
    n <- length(xs)
    X <- array(NA, c(n, q))
    X[, 1] <- 1
    X[, 2] <- xs
    X[, 3:q] <- outer(xs, xknot, FUN = basis)
    X
}

moveTo01 <- function(x) {
    x <- x - min(x)
    x <- x / max(x)
    x
}

makeXknot <- function(q) {
    xknot <- (1:(q-2))/(q-1)
    xknot
}

splineS <- function(xknot) {
    q <- length(xknot) + 2
    S <- array(0, c(q, q))
    S[3:q, 3:q] <- outer(xknot, xknot, FUN = basis)
    S
}

matRoot <- function(X) {
    d <- eigen(X, symmetric = TRUE)
    r <- d$vectors %*% diag(sqrt(d$values)) %*% t(d$vectors)
    r
}

penRegSpline <- function(y, xs, xknot, lambda) {
    q <- length(xknot) + 2
    n <- length(xs)
    Xa <- rbind(splineX(xs, xknot), matRoot(splineS(xknot)) * sqrt(lambda))
    y[(n+1):(n+q)] <- 0
    m <- lm(y ~ Xa - 1)
    m
}

plotPred <- function(y, xs, xknot, m) {
    xpred <- (0:100)/100
    Xp <- splineX(xpred, xknot)
    plot(xs, y, type='p')
    lines(xpred, Xp %*% coef(m))
}



## source('~/GitHub/pfLL/spline.R')
## load('~/GitHub/pfLL/pfLL.RData')

## xs <- moveTo01(out1$x)
## y <- out1$y

## q <- 6
## xknot <- makeXknot(q)
## X <- splineX(xs, xknot)
## X <- ns(xs, knots = xknot, intercept = TRUE)
## m <- lm(y ~ X - 1)
## plotPred(y, xs, xknot, m)

## xknot <- makeXknot(10)
## m2 <- penRegSpline(y, xs, xknot, lambda = 0.005)
## plotPred(y, xs, xknot, m2)







