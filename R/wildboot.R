##' Wild Bootstrap for regression model.
##'
##' Calculate several wild bootstrapped quantities.
##' 
##' @title wildboot
##' @param x an object
##' @param ... other arguments
##' @return A list with several components
##' @rdname wildboot
##' @author Giuseppe Ragusa
##' @export
wildboot <- function(x, ...)
    UseMethod('wildboot')

##' @method wildboot reg
##' @S3method wildboot reg
##' @return \code{NULL}
##' @rdname wildboot
##' @export
wildboot.reg <- function(obj, reps=999, null, 
                         type = c('radamacher', 'mtp', 'mn1', 'mn2'))
{

    ## null is a list
    ## for example reg(y~x)
    ## null = list(x=0)
    ## means we are testing whether
    ## the coefficient of x is 0
    ## if null is NULL
    ## then look whether
    ## t is a vector of dimension k
    ## t0 is numeric null hypothesis
    ## such that
    ## H_0: t%*%beta = t0
    ## In this case, the wildbootstrap does
    ## not impose the null
    types  <- c('radamacher', 'mtp', 'mn1', 'mn2')
    wbtype <- match.arg(type)
    wbtype <- switch(wbtype,
                     radamacher = 1,
                     mtp        = 2,
                     mn1        = 3,
                     mn2        = 4)
    
    z   <- obj
    y   <- z$model[[1]]
    X   <- model.matrix(z)
    b   <- coef(z)
    w   <- weights(z)
    n   <- length(y)
    k   <- ncol(X)
    r   <- residuals(z) 

    p <- z$rank    
    R <- chol2inv(z$qr$qr[1:p, 1:p, drop = FALSE])

    if(is.null(w))
        w <- rep(1, n)

    y <- y*sqrt(w)
    X <- X*c(sqrt(w))
    r <- r*sqrt(w)
    
    if(missing(null))
        stop("'null' is missing")

    if(!is.list(null))
        stop("'null' must be a list")
    ## Check null hypthesis        
    tmp    <- match(names(null), names(b), nomatch = NA)
    if(all(is.na(tmp))) {
        stop("'null' not defined properly")
    } else {
        wr     <- na.omit(tmp)[1]
        null   <- unlist(null)[1]   ## Null hypothesis
    }
    
    if(!is.null(z$cluster)) {
        cluster <- as.factor(z$cluster)
        j <- order(cluster)
        clus.size <- table(cluster)
        clus.start <- c(1, 1 + cumsum(clus.size))
        storage.mode(clus.start) <- "integer"
        nc <- length(levels(cluster))
        X <- X[j, , drop = FALSE]
        y <- y[j]
        clus.start <- clus.start[-(nc + 1)]        
    } else {
        nc <- 1
        clus.start <- c(1,n)
        clus.size <- c(n)
        fcl <- 1
    }    
    reps <- ceiling(reps)


    ## H0 Imposed

    ## Calculate beta_restricted and u_restricted
    yr  <- y-X[, wr, drop = FALSE]%*%null
    Xr  <- X[, -wr, drop = FALSE]
    br0 <- solve(crossprod(Xr), crossprod(Xr, yr)) ## this is ((k-1) x 1)
    Ur  <- yr-Xr%*%br0                              ## Restricted residuals
    br  <- rep(null, k)
    br[wr] <- null
    br[-wr] <- br0
    out <- .Call("wb_null", X, R, y, br, Ur, clus.start, clus.size, reps, wbtype, PACKAGE="grpack")
     
    coef.wb0            <- out$coef_wb[,wr]
    serr.wb0.HC1        <- out$sd_wb_HC1[,wr]
    serr.wb0.HC2        <- out$sd_wb_HC2[,wr]
    serr.wb0.HC3        <- out$sd_wb_HC3[,wr]

    names(coef.wb0) <- names(b)[wr]
    names(serr.wb0.HC1) <- names(b)[wr]
    names(serr.wb0.HC2) <- names(b)[wr]
    names(serr.wb0.HC3) <- names(b)[wr]

    tstat0.HC1 <- (coef.wb0 - null)/serr.wb0.HC1
    tstat0.HC2 <- (coef.wb0 - null)/serr.wb0.HC2
    tstat0.HC3 <- (coef.wb0 - null)/serr.wb0.HC3
    tstat0.se  <- (coef.wb0 - null)/sd(coef.wb0)
    
    ## Do the unconstrained wildbootstrap
    out <- .Call("wb_null", X, R, y, b, r, clus.start, clus.size, reps, wbtype, PACKAGE="grpack")
    coef.wb            <- out$coef_wb[,wr]
    serr.wb.HC1        <- out$sd_wb_HC1[,wr]
    serr.wb.HC2        <- out$sd_wb_HC2[,wr]
    serr.wb.HC3        <- out$sd_wb_HC3[,wr]
    
    names(coef.wb) <- names(b)[wr]
    names(serr.wb.HC1) <- names(b)[wr]
    names(serr.wb.HC2) <- names(b)[wr]
    names(serr.wb.HC3) <- names(b)[wr]

    tstat.HC1 <- (coef.wb - b[wr])/serr.wb.HC1
    tstat.HC2 <- (coef.wb - b[wr])/serr.wb.HC2
    tstat.HC3 <- (coef.wb - b[wr])/serr.wb.HC3    
    tstat.se  <- (coef.wb - b[wr])/sd(coef.wb0)

    tstat <- list(HC1 = tstat.HC1, HC2 = tstat.HC2, HC3 = tstat.HC3, se = tstat.se)
    tstat0 <- list(HC1 = tstat0.HC1, HC2 = tstat0.HC2, HC3 = tstat0.HC3, se = tstat0.se)

    obj$coef.wb       <- coef.wb        ## reps x 1
    obj$coef.wb0      <- coef.wb0       ## reps x 1

    obj$se.wb.HC1  <- serr.wb.HC1
    obj$se.wb.HC2  <- serr.wb.HC2
    obj$se.wb.HC3  <- serr.wb.HC3

    obj$se.wb0.HC1 <- serr.wb0.HC1
    obj$se.wb0.HC2 <- serr.wb0.HC2
    obj$se.wb0.HC3 <- serr.wb0.HC3

    obj$tstat0 <- tstat0
    obj$tstat  <- tstat
    
    obj$null         <- null         
    obj$reps         <- reps         
    obj$distr        <- types[wbtype]
    class(obj)       <- c("reg.wb", class(obj))
    obj$nc           <- nc
    obj$wr           <- wr
    obj
}

##' @method print reg.wb
##' @S3method print reg.wb
print.reg.wb <- function(x, ...) {
    summary.reg.wb(x)
    cat('Wild bootstrap:\n')
    cat('\n')
    cat(' # reps: ', x$reps, '; distr: ', x$distr, sep = '')
    if(!is.null(x$cluster))
        cat('; cluster:', paste(deparse(x$clusterby)))
    if(!is.null(x$weights))
        cat('; weighted:', paste(deparse(x$weightedby)))
    cat('\n\n')
        
    cat('Null hypotheses:\n')
    for(j in na.omit(x$wr))
        cat('', names(x$null), ' = ', x$null, '\n')
    cat('\n')
}

##' @method summary reg.wb
##' @S3method summary reg.wb
summary.reg.wb <- function(x, vcov. = "HC3", ...) {
    b  <- coef(x)
    se <- sqrt(diag(vcov(x, type = vcov.)))
    
    h0 <- unlist(x$null)
    tstat <- (b-h0)/se

    ats <- abs(tstat[x$wr])
    
    edf3 <- ecdf(x$tstat[[vcov.]])
    pv <- edf3(-abs(ats))+(1-edf3(abs(ats)))    
        
    edf3 <- ecdf(x$tstat0[[vcov.]])
    pv0 <- edf3(-abs(ats))+(1-edf3(abs(ats)))
    
    out <-  x
    out$pv0 <-  pv0
    out$pv  <-  pv
    out$vcov. <- vcov.
    class(out) <- c('summary.reg.wb', class(x))
    out
}

##' @method print summary.reg.wb
##' @S3method print summary.reg.wb
print.summary.reg.wb <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat('\nWild bootstrap: ')
    cat('# reps: ', x$reps, ', distr: ', x$distr, sep='')

    cat("\n\nModel:\n")
    cat(" ", paste(deparse(x$call), sep="\n", collapse = "\n"), "\n", sep="")
    if(!is.null(x$weights))
        cat('; weighted:', paste(deparse(x$weightedby)))
    cat('\n\n')
        
    cat('Null hypotheses:\n')
    cat('', names(x$null), ' = ', x$null, '\n')
    cat('\n')
    cat(' p-value (with Null imposed):', format.pval(x$pv0),
        symnum(x$pv0, corr = FALSE, na = FALSE, 
               cutpoints = c(0, 0.01, 0.05, 0.1, 1), 
               symbols = c("***", "**", "*", " ")), '\n')

    cat(' p-value (no null imposed):', format.pval(x$pv),
        symnum(x$pv, corr = FALSE, na = FALSE, 
               cutpoints = c(0, 0.01, 0.05, 0.1, 1), 
               symbols = c("***", "**", "*", " ")), '\n')

    cat("\n ---\n Signif. codes: '***' 0.01 '**' 0.05 '*' 0.1  \n")
    cat(" Variance type:", x$vcov., '\n')
}


## wildboot.reg2 <- function(object, sim = 999, null = list(),
##                           type = c('radamacher', 'mtp', 'mn1','mn2'),
##                           vcov = c('HC1', 'HC2', 'HC3'))
## {
##     type <- match.arg(type)
##     vcov <- match.arg(vcov)
##     ## Null hypothesis is a list
##     ## naming names
##     if(!is.list(null) || length(null)>1)
##         stop('null must be a list of max length 1')

##     mf <- model.frame(object)

##     y <- mf[,1]
##     X <- model.matrix(object)
##     b <- coef(object)
##     r <- residuals(object)
##     k <- length(b)
##     if(missing(null))
##         stop("'null' is missing")
    
##     if(!is.list(null))
##         stop("'null' must be a list")
##     ## Check null hypthesis        
##     tmp    <- match(names(null), names(b), nomatch = NA)
##     if(all(is.na(tmp))) {
##         stop("'null' not defined properly")
##     } else {
##         wr     <- na.omit(tmp)[1]
##         null   <- unlist(null)[1]   ## Null hypothesis
##     }
    
##     cluster <- as.factor(object$cluster)
##     j <- order(cluster)
##     clus.size <- table(cluster)
##     clus.start <- c(1, 1 + cumsum(clus.size))
##     nc <- length(levels(cluster))
##     clus.start <- clus.start[-(nc + 1)]
##     storage.mode(clus.start) <- "integer"
##     p <- object$rank
##     n <- NROW(object$qr$qr)
##     w <- object$weights
##     sp <- p
##     R <- chol2inv(object$qr$qr[1:p, 1:p, drop = FALSE])

##     factor <- sqrt((n-1)/(n-p) * nc/(nc-1))

##     y <- y[j]
##     X <- X[j, , drop = FALSE]
##     rr <- r[j]
    
##     out <- matrix(0, sim, 2*p)
##     ## weighting constants
##     tp1 <- -(sqrt(5)-1)/2
##     tp2 <- (sqrt(5)+1)/2
##     tpp <- (sqrt(5)+1)/(2*sqrt(5))
##     delta1 <- (3/4+sqrt(17)/12)^.5
##     delta2 <- (3/4-sqrt(17)/12)^.5
##     wbwf <- switch(type,
##                    radamacher = function(cs) {
##                        ww <- runif(nc)
##                        ww <- ifelse(ww<0.5,-1, 1)
##                        out <- NULL
##                        for(j in 1:length(ww)) out <- c(out, rep(ww[j], clus.size[j]))
##                        out
##                    },
##                    mtp = function(cs) {
##                        ww <- runif(nc)
##                        ww <- ifelse(ww<tpp,tp1, tp2)
##                        out <- NULL
##                        for(j in 1:length(ww)) out <- c(out, rep(ww[j], clus.size[j]))
##                        out
##                    },
##                    mn1 = function(cs) {
##                        ww <- rnorm(nc)
##                        ww <- ww/sqrt(2)+(ww^2-1)/2
##                        out <- NULL
##                        for(j in 1:length(ww)) out <- c(out, rep(ww[j], clus.size[j]))
##                        out
##                    },
##                    mn2 = function(cs) {
##                        ww <- (delta1+rnorm(nc)/sqrt(2))*(delta2+rnorm(nc)/sqrt(2))-delta1*delta2
##                        out <- NULL
##                        for(j in 1:length(ww)) out <- c(out, rep(ww[j], clus.size[j]))
##                        out
##                    })

##     mr2 <- function() {
##         res <- NULL
##         for (jj in 1:nc) {
##             ind   <- clus.start[jj]+ (0:(clus.size[jj]-1)) 
##             Xi    <- X[ind,,drop=FALSE]
##             Hgg   <- chol(diag(length(ind))-Xi%*%R%*%t(Xi), pivot = TRUE)
##             pivot <- attr(Hgg, "pivot")
##             oo    <- order(pivot)
##             Hgg   <- Hgg[,oo]
##             res   <- c(res, solve(Hgg)%*%r[ind])
##         }
##         res
##     }
    
##     mr3 <- function() {
##         res <- NULL
##         for (jj in 1:nc) {
##             ind <- clus.start[jj]+ (0:(clus.size[jj]-1)) 
##             Xi  <- X[ind,,drop=FALSE]
##             Hgg <- solve(diag(length(ind))-Xi%*%R%*% t(Xi))
##             res <- c(res, Hgg%*%r[ind])
##         }
##         sqrt((nc-1)/nc)*res
##     }

##     yr  <- y-X[, wr, drop = FALSE]%*%null
##     Xr  <- X[, -wr, drop = FALSE]
##     br0 <- solve(crossprod(Xr), crossprod(Xr, yr)) ## this is ((k-1) x 1)
##     Ur  <- yr-Xr%*%br0                              ## Restricted residuals
##     br  <- rep(null, k)
##     br[wr] <- null
##     br[-wr] <- br0
##     Yhat <- X%*%br
##     for(j in 1:sim)
##     {
##         wr <- -c(Ur)
##         yw <- Yhat+wr
##         betaw <- solve(crossprod(X), crossprod(X, yw))
##         r <- c(yw-X%*%betaw)
##         res <- switch(EXPR = vcov,                       
##                       HC2 = mr2(),
##                       HC3 = mr3(),
##                       factor*r
##                       )

##         score <- X*c(res)
##         W <- matrix(
##                     .Fortran("robcovf", n, sp, nc, clus.start, clus.size, 
##                              score, double(sp), double(sp * sp), w = double(sp * sp),
##                              PACKAGE = "grpack")$w, nrow = sp)
##         browser()
##         se <- sqrt(diag(R%*%W%*%t(R)))
##         out[j,] <- c(betaw, se)
##     }
    
##     ans <- list(coef.wb = out[,1:p], se.wb = out[,(p+1):(2*p)])
##     colnames(ans$boot.coef) <- rownames(betahat)
##     colnames(ans$ses) <- rownames(betahat)
    
##     ans$lm.full <- object
##     ans$lm.restricted <- rr
##     ans$coef <- betahat
##     ans$se <- sehat
##     attr(ans, 'clus.size') <- clus.size
##     attr(ans, 'clus.start') <- clus.start
##     attr(ans, 'which.restricted') <- whr
##     attr(ans,'restrictions') <- null
##     attr(ans,'p') <- p
##     class(ans) <- 'wildboot'
##     ans
## }



## ## wildbootreg <- function(obj, reps, null,
## ##                          type = c('radamacher', 'mtp', 'mn1', 'mn2'))
## ## {
## ##     wbtype <- match.arg(type)
## ##     wtyp <- switch(wbtype,
## ##                    radamacher = 1,
## ##                    mtp        = 2,
## ##                    mn1        = 3,
## ##                    mn2        = 4)
## ##     z   <- obj
## ##     y   <- z$model[[1]]
## ##     X   <- model.matrix(z)
## ##     b   <- coef(z)
## ##     n   <- NROW(y)
## ##     k   <- NCOL(X)
## ##     wr  <- 1:k
## ##     w   <- weights(z)

## ##     p <- z$rank    
## ##     R <- chol2inv(z$qr$qr[1:p, 1:p, drop = FALSE])
    
## ##     if(is.null(w))
## ##         w <- rep(1, n)
## ##     w <- n*w/sum(w)
## ##     res <- residuals(z)
    
## ##     if(missing(null) | !is.list(null))
## ##         stop('Null must be a list')
    
## ##     if(is.list(null)) {
## ##         tmp  <- match(names(b), names(null), nomatch = NA)
## ##         null <- unlist(null)[tmp]
## ##         wr   <- which(!is.na(tmp))
## ##     }
    
## ##     lwr <- sum(!is.na(wr))
    
## ##     if(missing(reps))
## ##         stop('reps must be an integer > 0')

## ##     reps <- ceiling(reps)
    
## ##     if(!is.null(z$cluster)) {
## ##         cluster <- as.factor(z$cluster)
## ##         j <- order(cluster)
## ##         clus.size <- table(cluster)
## ##         clus.start <- c(1, 1 + cumsum(clus.size))
## ##         nc <- length(levels(cluster))
## ##         X <- X[j, , drop = FALSE]
## ##         y <- y[j]
## ##         w <- w[j]
## ##         clus.start <- clus.start[-(nc + 1)]
## ##         fcl <- (nc/(nc-1))
## ##     } else {
## ##         nc <- 1
## ##         clus.start <- c(1,n)
## ##         clus.size <- c(n)
## ##         fcl <- 1
## ##     }
    
## ##     if(!is.null(w)) {
## ##         X      <- X*c(sqrt(w))
## ##         y      <- y*c(sqrt(w)) 
## ##     }

## ##     factor <- sqrt(((n-1)/(n-k))*fcl)
        
## ##     ## Containers
## ##     # These contains the t-test based on imposing the
## ##     # null hypothesis using a) serr.wb, and b) the
## ##     # se obtained by sd(coef.wb)
## ##     tstat0 <- tstat0.se  <- matrix(0, reps, lwr)
## ##     # These contains the t-test without imposing the
## ##     # null hypothesis using a) serr.wb and b) the
## ##     # seerr obtained by sd(coef.wb)
## ##     tstat  <- tstat.se    <- matrix(0, reps, k)
    
## ##     storage.mode(wtyp)       <- "integer"
## ##     storage.mode(clus.start) <- "integer"
## ##     storage.mode(clus.size)  <- "integer"
## ##     storage.mode(nc)         <- "integer"
## ##     storage.mode(n)          <- "integer"
## ##     storage.mode(k)          <- "integer"
## ##     storage.mode(reps)       <- "integer"
## ##     storage.mode(tmp)        <- "integer"
## ##     storage.mode(lwr)        <- "integer"
## ##     storage.mode(b)          <- "double"
## ##     storage.mode(factor)     <- "double"
## ##     storage.mode(X)          <- "double"
## ##     storage.mode(y)          <- "double"
## ##     storage.mode(res)        <- "double"
    
## ##     for(j in 1:lwr) {
## ##         tmp <- wr[j]
## ##         yr  <- y-X[, tmp, drop = FALSE]%*%null[[tmp]]
## ##         Xr  <- X[, -tmp, drop = FALSE]
## ##         br  <- solve(crossprod(Xr), crossprod(Xr, yr))
## ##         ur  <- yr-Xr%*%br 
## ##         kr  <- ncol(Xr)
## ##         brs <- rep(0,k)
## ##         brs[tmp] <- null[tmp]
## ##         brs[-tmp] <- br

## ##         storage.mode(brs) <- "double"
## ##         storage.mode(ur)  <- "double"
## ##         storage.mode(ur)  <- "double"
        
## ##         out <- .C('wildbootr', X, y, ur, brs, factor, tmp, lwr, n, k,
## ##                   reps, clus.start, clus.size, nc, wtyp,
## ##                   coef.wb = double(reps), serr.wb = double(reps))
        
## ##         coef.wb0        <- out$coef.wb
## ##         serr.wb0        <- out$serr.wb
## ##         names(coef.wb0) <- names(b)[tmp]
## ##         tstat0[,j]      <- (coef.wb0 - null[tmp])/serr.wb0
## ##         tstat0.se[,j]   <- (coef.wb0 - null[tmp])/sd(coef.wb0)
## ##     }
    
## ##     nul <- unlist(null)
## ##     nul[is.na(nul)] <- 0
## ##     tstat0.se <- (b[wr]-nul[wr])/sd(coef.wb0)
## ##     ## Do the unconstrained

## ##     out <- .C('wildbootr', X, y, res, b, factor, as.integer(0), k,
## ##               n, k, reps, clus.start, clus.size, nc, wtyp,
## ##               coef.wb = double(reps*k), serr.wb = double(reps*k))

## ##     coef.wb <- matrix(out$coef.wb, nrow = reps, ncol = k)
## ##     serr.wb <- matrix(out$serr.wb, nrow = reps, ncol = k)
## ##     colnames(coef.wb) <- names(b)
## ##     tms <- sd(coef.wb)    
## ##     for(j in 1:k) {
## ##         tmp <- (coef.wb[,j] - b[j])
## ##         tstat[,j]    <- tmp/serr.wb[,j]
## ##     }

## ##     tstat.se         <- (b-nul)/tms 
## ##     obj$coef.wb      <- coef.wb        ## reps x k
## ##     obj$se.wb        <- serr.wb        ## reps x k
## ##     obj$coef.wb0     <- coef.wb0       ## reps x lwr
## ##     obj$se.wb0       <- serr.wb0       ## reps x lwr
## ##     obj$tstat        <- tstat          ## reps x k
## ##     obj$tstat0       <- tstat0         ## reps x lwr
## ##     obj$tstat.se     <- tstat.se       ## k x 1
## ##     obj$tstat0.se    <- tstat0.se      ## lwr x 1
## ##     obj$null         <- null           ## lwr x 1
## ##     obj$wr           <- wr             ## k x 1 (NA non restricted)
## ##     obj$lwr          <- lwr
## ##     obj$se.wb.sd     <- tms            ## k x 1
## ##     obj$reps         <- reps           ## # of replications
## ##     obj$distr        <- 'em_a'         ## Only this right now
## ##     class(obj)       <- c("reg.wb", class(obj))
## ##     obj$model        <- NULL
## ##     obj$nc           <- nc
## ##     obj
## ## }
