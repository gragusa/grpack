##' \code{reg} is used to fit linear models. It can be used to carry
##' out regression, single stratum analysis of variance and analysis
##' of covariance (although aov may provide a more convenient
##' interface for these). It extends \code{lm} by allowing cluster
##' standard error and by defining a summary method which uses by
##' default heteroskedastic robust standard errors.
##'
##' Details are similar to \code{lm}.
##'
##' The only difference is that the factor at which level standard
##' errors are clustered must be specified.  
##' 
##' @title lm
##' @param formula an object of class \code{formula} (or one that can be
##' coerced to that class): a symbolic description of the model to be
##' fitted. See \code{lm} for details.
##' @param data an optional data frame, list or environment (or object
##' coercible by as.data.frame to a data frame) containing the
##' variables in the model. If not found in data, the variables are
##' taken from environment(formula), typically the environment from
##' which lm is called.
##' @param subset an optional vector specifying a subset of
##' observations to be used in the fitting process.
##' @param weights an optional vector of weights to be used in the
##' fitting process. Should be NULL or a numeric vector. If non-NULL,
##' weighted least squares is used with weights weights (that is,
##' minimizing sum(w*e^2)); otherwise ordinary least squares is
##' used. See \code{lm} for "details".
##' @param na.action a function which indicates what should happen
##' when the data contain NAs. The default is set by the na.action
##' setting of options, and is na.fail if that is unset. The
##' "factory-fresh" default is na.omit. Another possible value is
##' NULL, no action. Value na.exclude can be useful.
##' @param method the method to be used; for fitting, currently only
##' \code{method = "qr"} is supported; \code{method = "model.frame"}
##' returns the model frame (the same as with \code{model = TRUE}, see
##' below).
##' @param model  
##' @param x 
##' @param y 
##' @param qr 
##' @param singular.ok logical. If FALSE (the default in S but not in
##' R) a singular fit is an error.
##' @param contrasts an optional list. See the contrasts.arg of
##' model.matrix.default.
##' @param offset this can be used to specify an a priori known
##' component to be included in the linear predictor during
##' fitting. This should be NULL or a numeric vector of length equal
##' to the number of cases. One or more offset terms can be included
##' in the formula instead or as well, and if more than one are
##' specified their sum is used. See \code{model.offset}
##' @param cluster a factor specifying at which level clustering of
##' standard error should take place. 
##' @param ... additional arguments to be passed to the low level
##' regression fitting functions (see \code{lm}).
##' @return A list similar to the one returned by \code{lm}.
##' @author Giuseppe Ragusa
##' @export
reg <- function(formula, data, subset, weights, na.action,
                method = "qr", model = TRUE, x = FALSE, y = FALSE, qr = TRUE,
                singular.ok = TRUE, contrasts = NULL, offset, cluster, ...)
{     
    ret.x <- x
    ret.y <- y
    cl <- match.call()
    mf <- mfc <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
                 "offset", "cluster"), names(mf), 0L)
    ## m[7] gives cluster mfc[[m[7]]]
    ## m[4] gives weights
    ## m[3] gives subset
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (method == "model.frame") 
        return(mf)
    else if (method != "qr") 
        warning(gettextf("method = '%s' is not supported. Using 'qr'", 
            method), domain = NA)
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- as.vector(model.weights(mf))
    if (!is.null(w) && !is.numeric(w)) 
        stop("'weights' must be a numeric vector")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) == 1) 
            offset <- rep(offset, NROW(y))
        else if (length(offset) != NROW(y)) 
            stop(gettextf("number of offsets is %d, should equal %d (number of observations)", 
                length(offset), NROW(y)), domain = NA)
    }
    cluster <- model.cluster(mf)
    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = if (is.matrix(y)) matrix(, 0, 
            3) else numeric(0), residuals = y, fitted.values = 0 * 
            y, weights = w, rank = 0L, df.residual = if (is.matrix(y)) nrow(y) else length(y))
        if (!is.null(offset)) {
            z$fitted.values <- offset
            z$residuals <- y - offset
        }
    }
    else {
        x <- model.matrix(mt, mf, contrasts)
        z <- if (is.null(w)) 
            lm.fit(x, y, offset = offset, singular.ok = singular.ok, 
                ...)
        else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok, 
            ...)
    }
    class(z) <- c("reg", if (is.matrix(y)) "mlm", "lm")
    z$na.action <- attr(mf, "na.action")
    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- .getXlevels(mt, mf)
    z$call <- cl
    z$terms <- mt
    z$cluster   <- cluster
    z$clusterby <- mfc$cluster
    z$weightedby  <- mfc$weights
    z$subsetby  <- mfc$subset
    z$aliased   <- is.na(coef(z))
    if (model) 
        z$model <- mf
    if (ret.x) 
        z$x <- x
    if (ret.y) 
        z$y <- y
    if (!qr) 
        z$qr <- NULL
    z
}

##' @S3method summary reg
summary.reg <- function (object, type = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"),
                         correlation = FALSE, symbolic.cor = FALSE,
                         ...)
{            
    z <- object
    p <- z$rank
    type <- match.arg(type)
    
    if (p == 0) {
        r <- z$residuals
        n <- length(r)
        w <- z$weights
        if (is.null(w)) {
            rss <- sum(r^2)
        } else {
            rss <- sum(w * r^2)
            r <- sqrt(w) * r
        }
        resvar <- rss/(n - p)
        ans <- z[c("call", "terms")]
        class(ans) <- "summary.lm"
        ans$aliased <- is.na(coef(object))  # used in print method
        ans$residuals <- r
        ans$df <- c(0L, n, length(ans$aliased))
        ans$coefficients <- matrix(NA, 0L, 4L)
        dimnames(ans$coefficients)<-
            list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
        ans$sigma <- sqrt(resvar)
        ans$r.squared <- ans$adj.r.squared <- 0
        return(ans)
    }
    Qr <- object$qr
    if (is.null(z$terms) || is.null(Qr))
      stop("invalid \'lm\' object:  no 'terms' nor 'qr' component")
    n <- NROW(Qr$qr)
   
    if(!is.null(z$clusterby)) {
      rdf = length(unique(z$cluster))-1
    } else {
      rdf <- z$df.residual    
    }
    
    if(is.na(z$df.residual) || rdf != z$df.residual)
      warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    p1 <- 1:p
    ## do not want missing values substituted here
    r <- z$residuals
    f <- z$fitted.values
    w <- z$weights
    X <- model.matrix(object)
    if (is.null(w)) {
      mss <- if (attr(z$terms, "intercept"))
        sum((f - mean(f))^2) else sum(f^2)
      rss <- sum(r^2)
      w <- rep(1, n)
    } else {
      mss <- if (attr(z$terms, "intercept")) {
        m <- sum(w * f /sum(w))
        sum(w * (f - m)^2)
      } else sum(w * f^2)
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE]) ## solve(X'X)
    V <- vcov(object, type)
    se <- sqrt(diag(V[Qr$pivot[p1],Qr$pivot[p1]]))
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms")]
    ans$residuals <- r
    ans$coefficients <-
      cbind(est, se, tval, 2*pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(ans$coefficients)<-
      list(names(z$coefficients)[Qr$pivot[p1]],
           c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$aliased <- is.na(coef(object))  # used in print method
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
      df.int <- if (attr(z$terms, "intercept")) 1L else 0L
      ans$r.squared <- mss/(mss + rss)
      ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n - df.int)/rdf)
      if(type=='const')
        ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
                            numdf = p - df.int, dendf = rdf)
      else
        ans$fstatistic <- tryCatch(c(value = (coef(object)[-1]%*%solve(V[-1,-1])%*%coef(object)[-1])/(p-1),
                                     df = p-1), error = function(e) NULL)
      
    } else ans$r.squared <- ans$adj.r.squared <- 0
    ##    ans$cov.unscaled <- R
    ans$vcov <- V[Qr$pivot[p1],Qr$pivot[p1]]
    ans$se <- se
    dimnames(ans$vcov) <- dimnames(ans$coefficients)[c(1,1)]
    if (correlation) {
      ans$correlation <- (R * resvar)/outer(se, se)
      dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
      ans$symbolic.cor <- symbolic.cor
    }
    if(!is.null(z$na.action)) ans$na.action <- z$na.action
    ans$type <- type
    ans$clusterby <- z$clusterby
    ans$cluster <-  z$cluster
    ans$weightedby <- z$weightedby
    ans$weights <- z$weights
    class(ans) <- c("summary.reg","summary.lm")
    ans
}

##' @S3method print summary.reg
##' @export
print.summary.reg <-
    function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
	      signif.stars= getOption("show.signif.stars"),...)
{
    cat("\nModel:\n")#S: ' ' instead of '\n'
    cat(paste(deparse(x$call[[2]]), sep="\n", collapse = "\n"), "\n", sep="")
    if(!is.null(x$weightedby))
        cat("(weighted by", paste(deparse(x$weightedby)), " sum of wgt is ", format(sum(x$weights), digits=2), ")\n")
    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    ## cat(if(!is.null(x$w) && diff(range(x$w))) "Weighted ",
    ##     "Residuals:\n", sep="")
    ## if (rdf > 5L) {
    ##     nam <- c("Min", "1Q", "Median", "3Q", "Max")
    ##     rq <- if (length(dim(resid)) == 2)
    ##         structure(apply(t(resid), 1, quantile),
    ##     	      dimnames = list(nam, dimnames(resid)[[2]]))
    ##     else  structure(quantile(resid), names = nam)
    ##     print(rq, digits = digits, ...)
    ## }
    ## else if (rdf > 0L) {
    ##     print(resid, digits = digits, ...)
    ## } else { # rdf == 0 : perfect fit!
    ##     cat("ALL", df[1], "residuals are 0: no residual degrees of freedom!\n")
    ## }
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
        if (nsingular <- df[3L] - df[1L])
            cat("\nCoefficients: (", nsingular,
                " not defined because of singularities)\n", sep = "")
        else cat("\n")
        coefs <- x$coefficients
        if(!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames=list(cn, colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits=digits, signif.stars=signif.stars, na.print="NA", signif.legend = FALSE, ...)
        ##printMatCoef(coefs, digits=digits, signif.stars=signif.stars, na.print="NA", ...)
    }
    ##
    ##cat("\nResidual standard error:",
    ##	format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    if (!is.null(x$fstatistic)) {
        fpval <- if(length(x$fstatistic)==3)
            format.pval(pf(x$fstatistic[1L], x$fstatistic[2L],
                           x$fstatistic[3L], lower.tail = FALSE), digits=digits)
        else
            format.pval(pchisq(x$fstatistic[1L], x$fstatistic[2L],
                               lower.tail = FALSE), digits=digits)
        cat('---\n')
        cat("Multiple R-squared:", formatC(x$r.squared, digits=digits))
        cat(", Adjusted R-squared:",formatC(x$adj.r.squared,digits=digits),
            "\nF-statistic:", formatC(x$fstatistic[1], digits=digits),
            "on", x$fstatistic[2])
        if(x$type=='const') cat(" and", x$fstatistic[3])
        cat(" DF, p-value:", fpval, "\n")
    }
    if(x$type!="const") {
         if(is.null(x$clusterby))
             cat("Heteroskedastic robust std. err., type: ", x$type, '\n')
         else
             cat("Std err. adjusted for", length(unique(x$cluster))," clusters in", x$clusterby, "with", x$type, '\n')
     }
    cat("Signif. codes: ","'***' .001 '**' .01 '*' .05 '.' 0.1\n")
    correl <- x$correlation
    if (!is.null(correl)) {
	p <- NCOL(correl)
	if (p > 1L) {
	    cat("\nCorrelation of Coefficients:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
		print(symnum(correl, abbr.colnames = NULL))
	    } else {
                correl <- format(round(correl, 2), nsmall = 2, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
	}
    }
    cat("\n")#- not in S
    invisible(x)
}

##' @S3method vcov reg
vcov.reg <- function (object,
                      type = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"),
                      ...) {
    type <- match.arg(type)    
    iscluster <- !is.null(object$cluster)
    if(type=="const")
        iscluster <- FALSE
    z <- object
    Qr <- z$qr    
    if (is.null(z$terms) || is.null(Qr))
        stop("invalid \'lm\' object:  no 'terms' nor 'qr' component")
    p    <- z$rank
    p1   <- 1:p
    n    <- NROW(Qr$qr)
    rdf  <- n - p
    ##ind  <- !is.na(coef(z))
    nNA  <- Qr$pivot[p1]
    V    <- matrix(NA, length(Qr$pivot), length(Qr$pivot))
    R    <- chol2inv(Qr$qr[p1, p1, drop = FALSE]) ## solve(X'X)
    w    <- if(is.null(z$weights)) rep(1,n) else z$weights
    r    <- z$residuals
    if(!iscluster) {
        if(type=="const") {
                rss  <- sum(w * r^2)
                Vout <- R*rss/rdf
        } else {
            if(type!="HAC") {
                meat.  <- meatHC(object, type = type, ...)            
                meat.  <- meat.[nNA, nNA]
                bread. <- bread(object, type = type, ...)
                Vout   <- (1/nrow(estfun(object)))*(bread. %*% meat. %*% bread.)
            }
            else {
                Vout <- vcovHAC(object,  ...)
            }
        }
    } else {
        facj   <- n/(n-p)
        X          <- model.matrix(z)
        cluster    <- as.factor(z$cluster)
        nc         <- length(levels(cluster))
        j          <- order(cluster)
        clus.size  <- table(cluster)
        clus.start <- c(1, 1 + cumsum(clus.size))
        w <- w[j]
        r <- r[j]    
        X <- X[j, nNA, drop = FALSE]*c(sqrt(w))        
        r <- r*sqrt(w)
        
        ## mr2 <- function() {
        ##     res <- NULL
        ##     for (jj in 1:nc) {
        ##         ind   <- clus.start[jj]+ (0:(clus.size[jj]-1)) 
        ##         Xi    <- X[ind,,drop=FALSE]
        ##         Hgg   <- chol(diag(length(ind))-Xi%*%R%*%t(Xi), pivot = TRUE)
        ##         pivot <- attr(Hgg, "pivot")
        ##         oo    <- order(pivot)
        ##         Hgg   <- Hgg[,oo]
        ##         res   <- c(res, solve(Hgg)%*%r[ind])
        ##     }
        ##     res
        ## }
        
        ## mr3 <- function() {
        ##     res <- NULL
        ##     for (jj in 1:nc) {
        ##         ind <- clus.start[jj]+ (0:(clus.size[jj]-1)) 
        ##         Xi  <- X[ind,,drop=FALSE]
        ##         Hgg <- solve(diag(length(ind))-Xi%*%R%*% t(Xi))
        ##         res <- c(res, Hgg%*%r[ind])
        ##     }
        ##     sqrt((nc-1)/nc)*res
        ## }
        
        res <- switch(EXPR = type,                       
                      HC2 = .Call("resHC2", X, r, R, clus.start, clus.size, PACKAGE = "grpack"),
                      HC3 = .Call("resHC3", X, r, R, clus.start, clus.size, PACKAGE = "grpack"),
                      sqrt((n-1)/(n-p) * nc/(nc-1))*r
                      )
        
        score <- X*c(res)
        clus.start <- clus.start[-(nc + 1)]
        storage.mode(clus.start) <- "integer"
        sp <- p
        W <- matrix(
                    .Fortran("robcovf", n, sp, nc, clus.start, clus.size, 
                             score, double(sp), double(sp * sp), w = double(sp * sp),
                             PACKAGE = "grpack")$w, nrow = sp)
        Vout <- R%*%W%*%t(R)
    }
    V[nNA, nNA] <- Vout
    colnames(V) <- rownames(V) <- names(coef(object))
    V
}


coeftestdefault <- function (x, vcov. = NULL, df = NULL, ...) 
{
    est <- coef(x)
    if (is.null(vcov.)) 
        se <- vcov(x)
    else {
        if (is.function(vcov.)) 
            se <- vcov.(x)
        else se <- vcov.
    }
    se <- sqrt(diag(se))
    if (!is.null(names(est)) && !is.null(names(se))) {
        anames <- names(est)[names(est) %in% names(se)]
        est <- est[anames]
        se <- se[anames]
    }
    tval <- as.vector(est)/se
    if (is.null(df)) {
      if(!is.null(x$clusterby)) {
        rdf = length(unique(x$cluster))-1
      } else {
        rdf <- x$df.residual    
      }
    }
    
    if (df > 0) {
      pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
      cnames <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
      mthd <-ifelse(is.finite(df), "t", "z")
    }
    else {
      stop('df must >0')
    }
    rval <- cbind(est, se, tval, pval)
    colnames(rval) <- cnames
    class(rval) <- "coeftest"
    attr(rval, "method") <- paste(mthd, "test of coefficients")
    attr(rval, "df") <- df
    return(rval)
}

##' coef is a method  for performing z and (quasi-)t
##' tests of estimated coefficients through \code{reg}
##'
##' Details
##' @title Testing Estimated Coefficients
##' @param x a \code{reg} object
##' @param vcov a covariance type
##' @param df the degrees of freedom to be used. If this is a finite
##' positive number a t test with df degrees of freedom is
##' performed. In all other cases, a z test (using a normal
##' approximation) is performed. If the \code{reg} has a
##' \code{cluster} component and \df} is \code{NULL} a t test with G-1
##' degrees of freedom is performed. 
##' @param ... other arguments
##' @rdname coeftest
##' @S3method coeftest reg
coeftest.reg <- function(x, vcov.=c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"), df = NULL) {
    vcov. <- match.arg(vcov.)
    b <- coef(x)
    cluster <- FALSE
    
    if(is.null(df) & !is.null(x$clusterby)) {
        df <- length(unique(x$cluster))-1
        cluster <- TRUE
    } else {
      if(is.null(df)){
        df <- x$df.residual
      }
    }   
    cov <- vcov(x, type = vcov.)
    rval <- coeftestdefault(x, vcov.=cov, df = df)                    
    attr(rval, "vcov") <- vcov
    attr(rval, "df") <- df
    attr(rval, "cluster") <- cluster
    class(rval) <- c("coeftest.reg", "coeftest")
    return(rval)
}

##' @S3method print coeftest.reg
print.coeftest.reg <- function(x, ...)
{
  mthd <- attr(x, "method")
  if(is.null(mthd)) mthd <- "Test of coefficients"
  if(is.finite(df <- attr(x, "df")))
      cat(paste("\n", mthd," with ",df, " df:\n\n", sep = ""))
  else
      cat(paste("\n", mthd,":\n\n", sep = ""))
  printCoefmat(x, ...)
  cat("\n")
  invisible(x)
}


##' @S3method confint reg
confint.reg <- function (object, parm, level = 0.95,
                         vcov. = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"), df=NULL, ...) 
{
    ## ... other argument passsed down to vcov.reg (for instance when type is HAC)
    type <- match.arg(vcov.)
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    if(is.null(df)) {
      if(!is.null(object$clusterby)) {
        rdf <- length(unique(object$cluster))-1
      } else {
        rdf <- object$df.residual
      }  
    }
  
    pct <- format.perc(a, 3)
    fac <- qt(a, df = rdf)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
                                               pct))
    ses <- sqrt(diag(vcov(object, type = type, ...)))[parm]
    ci[] <- cf[parm] + ses %o% fac
    ci
}





## geboot.reg <- function(object, sim = 999, wbweights = c('radamacher', 'exp', 'mn1','mn2'), ...)
## {
##     type <- match.arg(wbweights)
##     ## Null hypothesis is a list
##     ## naming names
##     ##     if(!is.list(null) || length(null)>1)
##     ##         stop('null must be a list of max length 1')
##     mf <- model.frame(object)
##     y <- mf[,1]
##     X <- model.matrix(object)
##     betahat <- coef(object)
##     sehat <- sqrt(diag(vcov(object)))
##     nv <- names(betahat)
##     r <- residuals(object)
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
##     if(!is.null(w))
##         factor <- sqrt((sum(w)-1)/(sum(w)-p) * nc/(nc-1))
##     else
##         factor <- sqrt((n-1)/(n-p) * nc/(nc-1))
##     X <- X[j, , drop = FALSE]
##     y <- y[j]
##     r <- r[j]
##     out <- matrix(0, sim, 2*p)
##     ## weighting constants
##     tp1 <- -(sqrt(5)-1)/2
##     tp2 <- (sqrt(5)+1)/2
##     tpp <- (sqrt(5)+1)/(2*sqrt(5))
##     delta1 <- (3/4+sqrt(17)/12)^.5
##     delta2 <- (3/4-sqrt(17)/12)^.5
##     wbwf <- switch(type,
##                    radamacher = function() {
##                        ww <- runif(nc)
##                        ww <- ifelse(ww<0.5,0, 2) 
##                        rep(ww + 1 - sum(ww)/length(ww), clus.size)
##                    },
##                    mtp = function() {
##                        ww <- runif(nc)
##                        ww <- 1+ifelse(ww<tpp,tp1, tp2) 
##                        rep(ww + 1 - sum(ww)/length(ww), clus.size)
##                    },
##                    exp = function() {
##                        ww <- rexp(nc)
##                        rep(ww + 1 - sum(ww)/length(ww), clus.size)
##                    },
##                    mn1 = function() {
##                        ww <- rnorm(nc)
##                        ww <- 1+ww/sqrt(2)+(ww^2-1)/2
##                        rep(ww + 1 - sum(ww)/length(ww), clus.size)
##                    },
##                    mn2 = function() {
##                        ww <- 1+(delta1+rnorm(nc)/sqrt(2))*(delta2+rnorm(nc)/sqrt(2))-delta1*delta2
##                        rep(ww + 1 - sum(ww)/length(ww), clus.size)
##                    })
    
##     for(j in 1:sim)
##     {
##         ww <- wbwf()
##         #yw <- Yhat+wr
##         Xw <- X*c(ww)
##         XwX <- crossprod(Xw,X)
##         XwXi <- solve(XwX)
##         betaw <- solve(XwX, crossprod(Xw, y))
##         res <- factor*c(y-X%*%betaw)
##         score <- Xw*c(res)
##         W <- matrix(
##                     .Fortran("robcovf", n, sp, nc, clus.start, clus.size, 
##                              score, double(sp), double(sp * sp), w = double(sp * sp),
##                              PACKAGE = "grpack")$w, nrow = sp)
##         #se <- sqrt(diag(R%*%W%*%t(R)))
##         se <- sqrt(diag(XwXi%*%W%*%XwXi))
##         out[j,] <- c(betaw, se)
##     }
##     ans <- list(boot.coef = out[,1:p], ses = out[,(p+1):(2*p)])
##     colnames(ans$boot.coef) <- rownames(betahat)
##     colnames(ans$ses) <- rownames(betahat)
    
##     ans$lm.full <- object
##     ##    ans$lm.restricted <- rr
##     ans$coef <- betahat
##     ans$se <- sehat
##     attr(ans, 'clus.size') <- clus.size
##     attr(ans, 'clus.start') <- clus.start
##     ##    attr(ans, 'which.restricted') <- whr
##     ##    attr(ans,'restrictions') <- null
##     attr(ans,'p') <- p

##     class(ans) <- c('geboot', 'wildboot')
##     ans
## }

## plot.wildboot <- function(z, short = TRUE, trim = FALSE, qtrim = 0.01, ...)
## {
##     az <- attributes(z)
##     wr <- az$which.restricted
##     restrictions <- TRUE
##     p <- az$p
##     pp <- ceiling(az$p/4)
##     lmf <- z$lm.full
##     slmf <- summary(lmf)
##     if(length(az$restrictions)==0)
##         short <- FALSE
##     if(short)
##     {
##         bb <- z$boot.coef[,wr]
##         if(trim)
##             bbp <- trim(bb, q = qtrim)
##         else
##             bbp <- bb
##         tst <- (slmf$coef[wr,1]-az$restrictions[[1]])/slmf$coef[wr,2]
##         par(mfcol = c(1,2))
##         rb <- range(bbp)
##         cb <- coef(lmf)[wr]
##         sb <- coef(slmf)[wr,2]
##         rbb <- c(cb-4*sb,cb+4*sb)
##         lim <- c(min(rbb[1], rb[1]), max(rbb[2], rb[2])) 
##         xx <- seq(lim[1],lim[2], len = 100)
##         yy <- dnorm(xx,mean = cb, sd = sb)
        
##         hist(bbp, xlab = names(coef(z$lm.full))[wr],
##              main = 'Wild Bootstrap Coefficient',
##              probability = TRUE, ylim = c(0, max(yy)), xlim = lim)
##         ## Normal Approximation
##         lines(x = xx, y = yy, col = 'red', lty = 2)
##         tstat <- (bb-az$restrictions[[1]])/z$ses[,wr]
##         hist(tstat, xlab = names(coef(z$lm.full))[wr],
##              main = 'Wild Bootstrap t-statatistics',
##              probability = TRUE)
##         tmp <- quantile(tstat, p = c(0.025, 0.975))
##         points(y = 0, x = tmp[1], col = 'blue')
##         points(y = 0, x = tmp[2], col = 'blue')
##         abline(v = tmp[1], col = 'blue', lty = 2)
##         abline(v = tmp[2], col = 'blue', lty = 2)
##         points(y = 0, x = tst, col = 'red', lwd = 2)
##     }

## }


## summary.wildboot <- function(z, null = list(), ...)
## {
##     az <- attributes(z)
##     if(length(az$restrictions)>0)
##     {
##         wr <- az$which.restricted
##         btstat <- (z$boot.coef[,wr]-az$restrictions[[1]])/z$ses[,wr]
##         tstat <- (z$coef[wr]-az$restrictions[[1]])/z$se[wr]
##         edf <- ecdf(btstat)
##         pv <- edf(-abs(tstat))+(1-edf(abs(tstat)))
##         npv <- 2*(1-pnorm(abs(tstat)))
##         null <- az$restrictions
##         btstat <- (z$boot.coef[,wr]-z$coef[wr])/sd(z$boot.coef[,wr])
##         edf <- ecdf(btstat)
##         pv1 <- edf(-abs(tstat))+(1-edf(abs(tstat)))
##     }
##     else
##     {
##         nv <- names(z$beta)
##         if(!is.list(null) || lengh(null)>1)
##             stop('The model is unrestricted. "null" must be a list of max length 1')
##         for(j in 1:length(nv))
##             if(nv[j]==names(null))
##                 wr <- j
##                 btstat <- (z$boot.coef[,wr]-null[[1]])/z$ses[,wr]
##         tstat <- (z$coef[wr]-null[[1]])/z$se[wr]
##         edf <- ecdf(btstat)
##         pv <- edf(-abs(tstat))+(1-edf(abs(tstat)))
##         btstat <- (z$boot.coef[,wr]-z$coef[wr])/sd(z$boot.coef[,wr])
##         edf <- ecdf(btstat)
##         pv1 <- edf(-abs(tstat))+(1-edf(abs(tstat)))
##         npv <- 2*(1-pnorm(abs(tstat)))
##     }

##     ans <- list( lm.full = z$lm.full, restrictions = null,
##                 pvalue = pv, pvalue.1 = pv1,  normal.pvalue = npv)
##     class(ans) <- 'summary.wildboot'
##     ans
## }

## summary.geboot <- function(z, null, ...)
## {
##     az <- attributes(z)
##     nv <- names(z$coef)
##     if(missing(null))
##         stop('Null hyppthesis must be specified')
##     if(!is.list(null) || length(null)>1)
##         stop('The model is unrestricted. "null" must be a list of max length 1')
##     for(j in 1:length(nv))
##         if(nv[j]==names(null))
##             wr <- j
##     ses <- sqrt((z$boot.coef-z$coef)^2)
##     btstat1 <- ((z$boot.coef[,wr]-z$coef[wr])/z$ses[,wr])
##     btstat2 <- (z$boot.coef[,wr]-z$coef[wr])/z$se[wr]
##     tstat <- (z$coef[wr]-null[[1]])/z$se[wr]
##     edf <- ecdf(btstat1)
##     pv <- edf(-abs(tstat))+(1-edf(abs(tstat)))
##     edf <- ecdf(btstat2)
##     pv1 <- edf(-abs(tstat))+(1-edf(abs(tstat)))
##     npv <- 2*(1-pnorm(abs(tstat)))
##     ans <- list( lm.full = z$lm.full, restrictions = null,
##             pvalue = pv, pvalue.1 = pv1, normal.pvalue = npv)
##     class(ans) <- 'summary.geboot'
##     ans
## }

## print.summary.wildboot <- function(z, ...)
## {
##     cat('Unrestricted OLS\n')
##     lmf <- summary(z$lm.full)
##     print(lmf)
        
##     cat('\nNull Hypothesis: ')
##     cat(' ', names(z$restrictions), ' = ', z$restrictions[[1]], '\n\n')
##     cat(' Wild Boot  p-value:' , format.pval(z$pvalue),
##         symnum(z$pvalue, corr = FALSE, na = FALSE, 
##                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
##                symbols = c("***", "**", "*", ".", " ")), '\n')

##     cat(' Asy Normal p-value:' , format.pval(z$normal.pvalue),
##         symnum(z$normal.pvalue, corr = FALSE, na = FALSE, 
##                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
##                symbols = c("***", "**", "*", ".", " ")), '\n\n')
## }
    

## print.summary.geboot <- function(z, ...)
## {
##     cat('Unrestricted OLS\n')
##     lmf <- summary(z$lm.full)
##     print(lmf)
        
##     cat('\nNull Hypothesis: ')
##     cat(' ', names(z$restrictions), ' = ', z$restrictions[[1]], '\n\n')
##     cat(' Wild Boot  p-value:' , format.pval(z$pvalue),
##         symnum(z$pvalue, corr = FALSE, na = FALSE, 
##                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
##                symbols = c("***", "**", "*", ".", " ")), '\n')

##     cat(' Asy Normal p-value:' , format.pval(z$normal.pvalue),
##         symnum(z$normal.pvalue, corr = FALSE, na = FALSE, 
##                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
##                symbols = c("***", "**", "*", ".", " ")), '\n\n')
## }
    
