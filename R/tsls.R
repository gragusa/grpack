##################################################################
## tsls.R /
## Author: Giuseppe Ragusa 
## Time-stamp: "2011-11-28 17:37:38 gragusa" 
##
## Description:
##################################################################

ivreg <- function(y, ...)
  UseMethod("ivreg")


tsls2 <- function (y, X, Z, names=NULL, weights,
                   cluster=NULL, ...) {
  n <- length(y)
  p <- ncol(X)
  if(missing(weights)||is.null(weights)) {
    rep.weights <- FALSE
    weights <- rep(1,n) }
  else {
      if(!is.null(weights) && !is.numeric(weights) ) 
        stop("'weights' must be a numeric vector")
      if(length(weights)!=n)
        stop("'weights' not of right length")
      rep.weights <- TRUE
  }
  w <- weights/sum(weights)*n
  D <- diag(w)
  invZtZ <- solve(crossprod(D%*%Z,Z))
  XtZ <- crossprod(D%*%X, Z)
  V <- solve(XtZ %*% invZtZ %*% t(XtZ))
  A <- V %*% XtZ %*% invZtZ
  b <- A %*% crossprod(D%*%Z, y)
  residuals <- (y - X %*% b)
  ##   if(is.null(cluster))
  ##     V <- A %*% crossprod(D%*%Z*c(residuals)) %*% t(A)
  ##   else {
  ##     cluster <- as.factor(cluster)
  ##     j <- order(cluster)
  ##     score <- Z*c(residuals)*c(w)
  ##     score <- score[j, , drop = FALSE]
  ##     clus.size <- table(cluster)
  ##     clus.start <- c(1, 1 + cumsum(clus.size))
  ##     nc <- length(levels(cluster))
  ##     clus.start <- clus.start[-(nc + 1)]
  ##     storage.mode(clus.start) <- "integer"
  ##     sp <- NCOL(score)
  ##     W <- matrix(if (TRUE) 
  ##         .Fortran("robcovf", n, sp, nc, clus.start, clus.size, 
  ##             score, double(sp), double(sp * sp), w = double(sp * sp),
  ##                  PACKAGE = "grpack")$w
  ##     else .Fortran("robcovf", n, p, nc, clus.start, clus.size, 
  ##                   X, double(p), double(p * p), w = double(p * p))$w, nrow = sp)
  ##     V <- A%*%W%*%t(A)
  ##     clus.struct <- list()
  ##     clus.struct$start <- clus.start
  ##     clus.struct$size  <- clus.size
  ##     clus.struct$nclus  <- nc
  ##     attr(V,'cluster') <- clus.struct
  ##   }
  result<-list()
  result$n <- n
  result$p <- p
  b <- as.vector(b)
  names(b) <- names
  result$coefficients <- b
  rownames(V) <- colnames(V) <- names
  result$vcov<- V
  result$residuals <- as.vector(residuals)
  result$response <- y
  result$model.matrix <- X 
  result$instruments <- Z
  result$cluster <- cluster
  result$s2 <- sum(residuals^2)/(n-p)
  if(rep.weights)
      result$weights <- weights
  result
}

ivregdefault <- function(y, X, Z, names=colnames(X), weights,
                         method = "tsls", cluster=cluster, start, ...)
{
    if(missing(weights))
        weights <- rep(1, length(y))
    if(any(weights)<0)
        stop('negative weights are not allowed')
    
    if(missing(cluster))
        cluster <- rep(1, length(y))
    
    if((length(y)!=NROW(X))|(length(y)!=NROW(Z))|(NROW(Z)!=NROW(X)))
        stop('mismetch')
    
    if(NCOL(X)>NCOL(Z))
        stop('order condition is not satisfied, NCOL(X)>NCOL(Z)')
    
    if(!is.function(method)&!is.character(method))
        stop("'method' is wrong")

    args <- list(...)
    
    lm <- try(lowerize(method), silent = TRUE)
    
    if(is.null(start) & lm != 'tsls')
        start <- tsls2(y,X,Z,names=names, weights=weights, cluster = cluster)$coef
    
    momfiv <- function(b) 
        Z*c(y-X%*%b)
    dmomfiv <- function(b, weights = 1) 
        crossprod(Z * c(weights), -X)
    data <- list(y = y, X = X, Z = Z)
    
    out <- if(is.function(method) && !is.null(method()$family)) 
        mdest(momfiv, start = start,
              dmomfun = dmomfiv, weights = weights, mdfamily = method, ...)
    else 
        switch(lowerize(method),
               tsls = tsls2(y,X,Z,names = names, weights = weights,
               cluster = cluster),
               liml = tsls2(y,X,Z,names = names, weights = weights,
               cluster = cluster),
               gmm = gmmest(momfiv, start = start, data = data,
                 dmomfun = dmomfiv, weights = weights, ...),
               el = mdest(momfiv, start = start, data = data,
                 dmomfun = dmomfiv, weights = weights, mdfamily = mdel, ...),
               et = mdest(momfiv, start = start, data = data,
                 dmomfun = dmomfiv, weights = weights, mdfamily = mdet, ...),
               cue = mdest(momfiv, start = start, data = data,
                 dmomfun = dmomfiv, weights = weights, mdfamily = mdcue, ...),
               ht = mdest(momfiv, start = start, data = data,
                 dmomfun = dmomfiv, weights = weights, mdfamily = mdht, ...),
               qt = mdest(momfiv, start = start, data = data,
                 dmomfun = dmomfiv, weights = weights, mdfamily = mdqt, ...)
             )
}    


ivreg <- function (formula, instruments, data, subset, weights,
                   na.action, contrasts = NULL,
                   cluster = NULL, method = "tsls", start = NULL,...) {
    if (missing(na.action))
        na.action <- options()$na.action
    m <- match.call(expand.dots = FALSE)
    
    mf <- match(c("formula", "instruments", "data", "subset", "weights", "na.action", 
                  "contrast", "cluster"), names(m), 0L)
    m <- m[c(1L, mf)]
    m$drop.unused.levels <- TRUE
    if (is.matrix(eval(m$data, sys.frame(sys.parent()))))
        m$data <- as.data.frame(data)
    response.name <- deparse(formula[[2]])
    form <- as.formula(paste(
                             paste(response.name, collapse = ""),
                             "~",
                             paste(deparse(formula[[3]]), collapse = ""),
                             "+",
                             paste(deparse(instruments[[2]]), collapse = "")))
    m$formula <- form
    m$instruments <- m$contrasts <- NULL
    m$method <- NULL
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, sys.frame(sys.parent()))
    na.act <- attr(mf, "na.action")
    Z <- model.matrix(instruments, data = mf, contrasts)
    y <- mf[, response.name]
    X <- model.matrix(formula, data = mf, contrasts)
    weights <- model.weights(mf)
    cluster <- model.cluster(mf)
    result <- ivregdefault(y, X, Z, names = colnames(X), weights = weights,
                           method = method, cluster = cluster, start, ...)
    result$response.name <- response.name
    result$formula <- formula
    result$instruments <- Z
    result$formula.Z <- instruments
    result$call <- match.call(expand.dots = TRUE)
    result$terms <- attr(mf, "terms")
    if (!is.null(na.act))
        result$na.action <- na.act
    if(!is.null(result$cluster))
        attr(result$vcov,'cluster')$name <- as.character(match.call()['cluster'])
    class(result) <- c(class(result), 'ivreg', 'grpack')
    result
}


print.ivreg <- function(x, ...) {
    cat("\nModel Formula: ")
    print(x$formula)
    cat("\nInstruments: ")
    print(x$formula.Z)
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\n")
    invisible(x)
}

summary.ivreg <- function(object, digits = 4, ...) {
    save.digits <- unlist(options(digits = digits))
    on.exit(options(digits = save.digits))
    cat("\n 2SLS Estimates (Robust Var)\n")
    if(!is.null(AC <- attr(object$vcov,'cluster')))
      {
          cat("\nNumber of cluster (", AC$name, "): ", AC$nclus,"\n")
      }
    cat("\nModel Formula: ")
    print(object$formula)
    cat("\nInstruments: ")
    print(object$formula.Z)
    cat("\nResiduals:\n")
    print(summary(residuals(object)))
    cat("\n")
    df <- object$n - object$p
    std.errors <- sqrt(diag(object$vcov))
    b <- object$coefficients
    t <- b/std.errors
    p <- 2*(1 - pt(abs(t), df))
    table <- cbind(b, std.errors, t, p)
    rownames(table) <- names(b)
    colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)")
    print(table)
    cat(paste("\nResidual standard error:", round(object$s, digits),
              "on", df, "degrees of freedom\n\n"))
}

## wildboot

disp <- function(x, ...)
    UseMethod('disp')

disp.default <- function(x, alpha = 0.05, ...)
{
    out <- matrix(0, p <- NROW(coef(x)), 2)
    out[,1] <- coef(x)
    out[,2] <- se(x)
    colnames(out) <- c('Estimate', 'Std. Error')
    rownames(out) <- names(coef(x))
    t <- abs(coef(x)/se(x))
    ind <- (t > abs(qnorm(alpha/2)))
    out <- out[ind, , drop = FALSE]
    rownames(out) <- names(coef(x))[ind]
    cat("\n")
    print(out, print.gap = 2)
    cat("\n")
}


se <- function(x, ...)
    UseMethod('se')

se.ivreg <- function(x)
    sqrt(diag(x$vcov))

se.reg <- function(x)
    sqrt(diag(vcov(x)))

se.md <- function(x, vcov. = 'wr')
    sqrt(diag(x$V[[vcov.]]))

model.cluster <- function (x) 
  x$"(cluster)"




vcov.ivreg <- function (object, hc = c('HC', 'HC1', 'HC2', 'HC3', 'iid'), ...)
{
    z <- object
    hc <- match.arg(hc)
    iscluster <- !is.null(z$cluster)
    if(iscluster & hc == 'HC1')
        hc = 'HC'
    ## do not want missing values substituted here
    r <- z$residuals
    #f <- z$fitted.values
    w <- z$weights
    X <- z$model.matrix
    Z <- z$instruments
    k <- NCOL(X)
    p <- NCOL(Z)
    n <- z$n
    if (is.null(w))
        w <- rep(1, n)
    else 
        r <- sqrt(w) * r
    XZ <- crossprod(X,Z)
    ZZ <- crossprod(Z)
    A <- solve(XZ%*%solve(ZZ, t(XZ)))
    ZZi <- solve(ZZ)
    R <- A%*%XZ%*%ZZi
    
    if(iscluster)
    {
        cluster <- as.factor(z$cluster)
        j <- order(cluster)
        clus.size <- table(cluster)
        clus.start <- c(1, 1 + cumsum(clus.size))
        nc <- length(levels(cluster))
        X <- X[j, , drop = FALSE]
        Z <- Z[j, , drop = FALSE]
        r <- r[j]
        mr2 <- function()
        {
            res <- NULL
            for (jj in 1:nc) {
                ind <- clus.start[jj]:(clus.start[jj+1]-1)
                Hgg <- chol(diag(length(ind))-Z[ind,,drop=FALSE]%*%ZZi%*% t(Z[ind,,drop=FALSE]))
                res <- c(res, solve(Hgg)%*%r[ind])
            }
            res
        }

        mr3 <- function()
        {
            res <- NULL
            for (jj in 1:nc) {
                ind <- clus.start[jj]:(clus.start[jj+1]-1)
                Hgg <- solve(diag(length(ind))-Z[ind,,drop=FALSE]%*%ZZi%*% t(Z[ind,,drop=FALSE]))
                res <- c(res, solve(Hgg)%*%r[ind])
            }
            res
        }
        
                
        res <- switch(hc,
                      HC  = sqrt((sum(w)-1)/(sum(w)-k) * nc/(nc-1))*r,
                      HC2 = mr2(),
                      HC3 = sqrt((nc-1)/nc)*mr3())
        
        score <- Z*c(res)
        clus.start <- clus.start[-(nc + 1)]
        storage.mode(clus.start) <- "integer"
        sp <- p
        W <- matrix(
                    .Fortran("robcovf", n, sp, nc, clus.start, clus.size, 
                             score, double(sp), double(sp * sp), w = double(sp * sp),
                             PACKAGE = "grpack")$w, nrow = sp)
        V <- R%*%W%*%t(R)
    }
    
    if(!iscluster)
        V <- switch(hc,
                    iid = {
                        rss <- sum(r^2)
                        resvar <- rss/(n-p)
                        V <- A * resvar
                        se <- sqrt(diag(V))},
                    HC = R%*%crossprod(Z*r)%*%t(R),
                    HC1 = n/(n-k)*R%*%crossprod(Z*c(r))%*%t(R),
                    HC2 = {
                        ktt <- diag(Z%*%ZZi%*%t(Z))
                        nr <- r/(1-ktt)
                        n/(n-p)*R%*%crossprod(Z*c(r))%*%t(R)},
                    HC3 = {
                        ktt <- diag(Z%*%ZZi%*%t(Z))
                        us <- r/(1-ktt)
                        ((n-1)/n)*R%*%(t(Z)%*%diag(us^2)%*%Z-(1/n)*(crossprod(Z*c(us))))%*%t(R)
                    })

    colnames(V) <- rownames(V) <- names(coef(object))
    V
}
