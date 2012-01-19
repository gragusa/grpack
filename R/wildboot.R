##' @export
wildboot.reg <- function(obj, reps, null,
                         wbtype = c('radamacher', 'mtp', 'mn1', 'mn2'))
{
    wbtype <- match.arg(wbtype)
    wtyp <- switch(wbtype,
                   radamacher = 1,
                   mtp        = 2,
                   mn1        = 3,
                   mn2        = 4)
    z   <- obj
    y   <- z$model[[1]]
    X   <- model.matrix(z)
    b   <- coef(z)
    n   <- NROW(y)
    k   <- NCOL(X)
    wr  <- 1:k
    w   <- weights(z)
    if(is.null(w))
        w <- rep(1, n)
    w <- n*w/sum(w)
    res <- residuals(z)
    
    if(missing(null) | !is.list(null))
        stop('For now')
    
    if(is.list(null))
    {
        tmp  <- match(names(b), names(null), nomatch = NA)
        null <- unlist(null)[tmp]
        wr   <- which(!is.na(tmp))
    }
    
    lwr <- sum(!is.na(wr))
    
    if(missing(reps))
        stop('reps must be an integer > 0')

    reps <- ceiling(reps)
    
    if(!is.null(z$cluster)) {
        cluster <- as.factor(z$cluster)
        j <- order(cluster)
        clus.size <- table(cluster)
        clus.start <- c(1, 1 + cumsum(clus.size))
        nc <- length(levels(cluster))
        X <- X[j, , drop = FALSE]
        y <- y[j]
        w <- w[j]
        clus.start <- clus.start[-(nc + 1)]
        fcl <- (nc/(nc-1))
    }
    else {
        nc <- 1
        clus.start <- c(1,n)
        clus.size <- c(n)
        fcl <- 1
    }
    
    if(!is.null(w))
    {
        X      <- X*c(sqrt(w))
        y      <- y*c(sqrt(w)) 
    }

    factor <- sqrt(((n-1)/(n-k))*fcl)
        
    ## Containers
    # These contains the t-test based on imposing the
    # null hypothesis using a) serr.wb, and b) the
    # se obtained by sd(coef.wb)
    tstat0 <- tstat0.se  <- matrix(0, reps, lwr)
    # These contains the t-test without imposing the
    # null hypothesis using a) serr.wb and b) the
    # seerr obtained by sd(coef.wb)
    tstat  <- tstat.se    <- matrix(0, reps, k)
    
    storage.mode(wtyp)       <- "integer"
    storage.mode(clus.start) <- "integer"
    storage.mode(clus.size)  <- "integer"
    storage.mode(nc)         <- "integer"
    storage.mode(n)          <- "integer"
    storage.mode(k)          <- "integer"
    storage.mode(reps)       <- "integer"
    storage.mode(tmp)        <- "integer"
    storage.mode(lwr)        <- "integer"
    storage.mode(b)          <- "double"
    storage.mode(factor)     <- "double"
    storage.mode(X)          <- "double"
    storage.mode(y)          <- "double"
    storage.mode(res)        <- "double"
    
    for(j in 1:lwr)
    {
        tmp <- wr[j]
        yr  <- y-X[, tmp, drop = FALSE]%*%null[[tmp]]
        Xr  <- X[, -tmp, drop = FALSE]
        br  <- solve(crossprod(Xr), crossprod(Xr, yr))
        ur  <- if(is.null(w)) yr-Xr%*%br else (yr-Xr%*%br)*c(w)
        kr  <- ncol(Xr)
        brs <- rep(0,k)
        brs[tmp] <- null[tmp]
        brs[-tmp] <- br

        storage.mode(brs) <- "double"
        storage.mode(ur) <- "double"
        storage.mode(ur)  <- "double"
        
        out <- .C('wildbootr', X, y, ur, brs, factor, tmp, lwr, n, k,
                  reps, clus.start, clus.size, nc, wtyp,
                  coef.wb = double(reps), serr.wb = double(reps))
        
        coef.wb0        <- out$coef.wb
        serr.wb0        <- out$serr.wb
        names(coef.wb0) <- names(b)[tmp]
        tstat0[,j]      <- (coef.wb0 - null[tmp])/serr.wb0
        tstat0.se[,j]   <- (coef.wb0 - null[tmp])/sd(coef.wb0)
    }
    nul <- unlist(null)
    nul[is.na(nul)] <- 0
    tstat0.se <- (b[wr]-nul[wr])/sd(coef.wb0)
    ## Do the unconstrained

    out <- .C('wildbootr', X, y, res, b, factor, as.integer(0), k,
              n, k, reps, clus.start, clus.size, nc, wtyp,
              coef.wb = double(reps*k), serr.wb = double(reps*k))

    coef.wb <- matrix(out$coef.wb, nrow = reps, ncol = k)
    serr.wb <- matrix(out$serr.wb, nrow = reps, ncol = k)
    colnames(coef.wb) <- names(b)
    tms <- sd(coef.wb)    
    for(j in 1:k) {
        tmp <- (coef.wb[,j] - b[j])
        tstat[,j]    <- tmp/serr.wb[,j]
    }

    tstat.se         <- (b-nul)/tms 
    obj$coef.wb      <- coef.wb        ## reps x k
    obj$se.wb        <- serr.wb        ## reps x k
    obj$coef.wb0     <- coef.wb0       ## reps x lwr
    obj$se.wb0       <- serr.wb0       ## reps x lwr
    obj$tstat        <- tstat          ## reps x k
    obj$tstat0       <- tstat0         ## reps x lwr
    obj$tstat.se     <- tstat.se       ## k x 1
    obj$tstat0.se    <- tstat0.se      ## lwr x 1
    obj$null         <- null           ## lwr x 1
    obj$wr           <- wr             ## k x 1 (NA non restricted)
    obj$lwr          <- lwr
    obj$se.wb.sd     <- tms            ## k x 1
    obj$reps         <- reps           ## # of replications
    obj$distr        <- 'em_a'         ## Only this right now
    class(obj)       <- c("reg.wb", class(obj))
    obj$model        <- NULL
    obj$nc           <- nc
    obj
}
##' @S3method print reg.wb
print.reg.wb <- function(x, ...) {
    print.reg(x)
    cat('Wild bootstrap:\n')
    cat('\n')
    cat(' # reps: ', x$reps, '; distr: ', x$distr, sep = '')
    if(!is.null(x$cluster))
        cat('; cluster:', x$clusterby)
    if(!is.null(x$weights))
        cat('; weighted:', x$weightedby)
    cat('\n\n')
        
    cat('Null hypotheses:\n')
    for(j in na.omit(x$wr))
        cat('', names(x$null)[j], ' = ', x$null[[j]], '\n')
    cat('\n')
}
##' @S3method summary reg.wb
summary.reg.wb <- function(x, ...) {
    b  <- coef(x)
    se <- sqrt(diag(vcov(x, ...)))
    h0 <- unlist(x$null)
    h0[is.na(h0)] <- 0
    tstat <- (b-h0)/se
    pv0 <- matrix(0, 2, x$lwr)
    pv  <- matrix(0, 2, length(b))

    for(j in 1:length(b))
    {   ats <- abs(tstat[j])
        edf3 <- ecdf(x$tstat[,j])
        pv[1,j] <- edf3(-abs(ats))+(1-edf3(abs(ats)))
    }
    pv[2,] <- 2*pnorm(-abs(x$tstat.se))
        
    for(j in 1:x$lwr)
    {   ts <- abs(tstat[x$wr[j]])
        edf3 <- ecdf(x$tstat0[,j])
        pv0[1,j] <- edf3(-ats)+(1-edf3(ats))
    }
    pv0[2,] <- 2*pnorm(-abs(x$tstat0.se))
    
    out              <-  x
    out$pv.tstat0.se <-  pv0[2,]
    out$pv.tstat.se  <-  pv[2,]
    out$pv.tstat     <-  pv[1,]
    out$pv.tstat0    <-  pv0[1,]
    if(!is.null(x$cluster) & x$nc>2) {
        out$pv.tstat.se.t  <- 2*pt(-abs(x$tstat.se),  df = x$nc-2)
        out$pv.tstat0.se.t <- 2*pt(-abs(x$tstat0.se), df = x$nc-2)
    }
    class(out) <- c('summary.reg.wb', class(x))
    out
}


##' @S3method print summary.reg.wb
print.summary.reg.wb <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    b <- coef(x)
    nb <- names(b)
    k <- length(b)
    cat('\nNull Hypothesis:')
    for(j in 1:k) {
        cat('\n ', nb[j], ' = ', x$null[j], '\n')
        cat(' wb p-value (a):' , format.pval(x$pv.tstat0[j]),
            symnum(x$pv.tstat0[j], corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                   symbols = c("***", "**", "*", ".", " ")), '\n')
        cat(' wb p-value (b):' , format.pval(x$pv.tstat0.se[j]),
            symnum(x$pv.tstat0.se[j], corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                   symbols = c("***", "**", "*", ".", " ")), '\n')
        cat(' wb p-value (c):' , format.pval(x$pv.tstat[j]),
            symnum(x$pv.tstat[j], corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                   symbols = c("***", "**", "*", ".", " ")), '\n')
        cat(' wb p-value (d):' , format.pval(x$pv.tstat.se[j]),
            symnum(x$pv.tstat.se[j], corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                   symbols = c("***", "**", "*", ".", " ")), '\n')
    }

    cat("\n ---\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
}



wildboot.reg2 <- function(object, sim = 999, null = list(),
                         wbweights = c('radamacher', 'mtp', 'mn1','mn2'), ...)
{
    type <- match.arg(wbweights)
    ## Null hypothesis is a list
    ## naming names
    if(!is.list(null) || length(null)>1)
        stop('null must be a list of max length 1')
    mf <- model.frame(object)
    y <- mf[,1]
    X <- model.matrix(object)
    betahat <- coef(object)
    sehat <- sqrt(diag(vcov(object)))
    nv <- names(betahat)
    rr <- object
    whr <- NA
    if(length(null)>0)
    {
        form <- paste(object$term[[2]], '~ 1')
        for(j in 1:length(nv))
        {
            if(nv[j]==names(null))
            {
                offset = null[[1]]*mf[[nv[j]]]
                whr <- j
            }
            else
                if(nv[j]!='(Intercept)')
                    form <- paste(form, nv[j], sep = '+') 
        }
        formula <- as.formula(form)
        rr$call[[2]] <- form
        rr$call$offset <- offset
        rr <-  eval(rr$call, envir = attributes(rr$terms)$.Environment)
    }
    se <- sqrt(diag(vcov(rr)))
    Yhat <- predict(rr)
    r <- residuals(rr)
    
    cluster <- as.factor(rr$cluster)
    j <- order(cluster)
    clus.size <- table(cluster)
    clus.start <- c(1, 1 + cumsum(clus.size))
    nc <- length(levels(cluster))
    clus.start <- clus.start[-(nc + 1)]
    storage.mode(clus.start) <- "integer"
    p <- object$rank
    n <- NROW(object$qr$qr)
    w <- object$weights
    sp <- p
    R <- chol2inv(object$qr$qr[1:p, 1:p, drop = FALSE])
    if(!is.null(w))
        factor <- sqrt((sum(w)-1)/(sum(w)-p) * nc/(nc-1))
    else
        factor <- sqrt((n-1)/(n-p) * nc/(nc-1))
    X <- X[j, , drop = FALSE]
    y <- y[j]
    r <- r[j]
    
    out <- matrix(0, sim, 2*p)
    ## weighting constants
    tp1 <- -(sqrt(5)-1)/2
    tp2 <- (sqrt(5)+1)/2
    tpp <- (sqrt(5)+1)/(2*sqrt(5))
    delta1 <- (3/4+sqrt(17)/12)^.5
    delta2 <- (3/4-sqrt(17)/12)^.5
    wbwf <- switch(type,
                   radamacher = function() {
                       ww <- runif(nc)
                       ww <- ifelse(ww<0.5,-1, 1) 
                       rep(ww, each = clus.size)
                   },
                   mtp = function() {
                       ww <- runif(nc)
                       ww <- ifelse(ww<tpp,tp1, tp2) 
                       rep(ww, each = clus.size)
                   },
                   mn1 = function() {
                       ww <- rnorm(nc)
                       ww <- ww/sqrt(2)+(ww^2-1)/2
                       rep(ww, each = clus.size)
                   },
                   mn2 = function() {
                       ww <- (delta1+rnorm(nc)/sqrt(2))*(delta2+rnorm(nc)/sqrt(2))-delta1*delta2
                       wr <- rep(ww, each = clus.size)
                   })
    
    for(j in 1:sim)
    {
        wr <- wbwf()*c(r)
        yw <- Yhat+wr
        betaw <- solve(crossprod(X), crossprod(X, yw))
        res <- factor*c(yw-X%*%betaw)
        score <- X*c(res)
        W <- matrix(
                    .Fortran("robcovf", n, sp, nc, clus.start, clus.size, 
                             score, double(sp), double(sp * sp), w = double(sp * sp),
                             PACKAGE = "grpack")$w, nrow = sp)
        se <- sqrt(diag(R%*%W%*%t(R)))
        out[j,] <- c(betaw, se)
    }
    ans <- list(coef.wb = out[,1:p], se.wb = out[,(p+1):(2*p)])
    colnames(ans$boot.coef) <- rownames(betahat)
    colnames(ans$ses) <- rownames(betahat)
    
    ans$lm.full <- object
    ans$lm.restricted <- rr
    ans$coef <- betahat
    ans$se <- sehat
    attr(ans, 'clus.size') <- clus.size
    attr(ans, 'clus.start') <- clus.start
    attr(ans, 'which.restricted') <- whr
    attr(ans,'restrictions') <- null
    attr(ans,'p') <- p
    class(ans) <- 'wildboot'
    ans
}
