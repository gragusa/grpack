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
    z$weighdby  <- mfc$weights
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


print.reg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    if(length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits=digits),
                      print.gap = 2, quote = FALSE)
    } else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

summary.reg <- function (object, correlation = FALSE, symbolic.cor = FALSE,
                         hc = c('HC', 'HC1', 'HC2', 'HC3', 'iid'), ...)
{
    hc <- match.arg(hc)
    z <- object
    p <- z$rank
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
    rdf <- n - p
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
    V <- vcov(object, hc)
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
        if(hc=='iid')
            ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
                                numdf = p - df.int, dendf = rdf)
        else
            ans$fstatistic <- tryCatch(c(value = coef(object)[-1]%*%solve(V[-1,-1])%*%coef(object)[-1],
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
    class(ans) <- c("summary.reg","summary.lm")
    ans
}

print.summary.reg <-
    function (x, digits = max(3, getOption("digits") - 3),
              symbolic.cor = x$symbolic.cor,
	      signif.stars= getOption("show.signif.stars"),...)
{
    cat("\nCall:\n")#S: ' ' instead of '\n'
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    resid <- x$residuals
    df <- x$df
    rdf <- df[2L]
    cat(if(!is.null(x$w) && diff(range(x$w))) "Weighted ",
        "Residuals:\n", sep="")
    if (rdf > 5L) {
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2)
	    structure(apply(t(resid), 1, quantile),
		      dimnames = list(nam, dimnames(resid)[[2]]))
	else  structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
    }
    else if (rdf > 0L) {
	print(resid, digits = digits, ...)
    } else { # rdf == 0 : perfect fit!
	cat("ALL", df[1], "residuals are 0: no residual degrees of freedom!\n")
    }
    if (length(x$aliased) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
        if (nsingular <- df[3L] - df[1L])
            cat("\nCoefficients: (", nsingular,
                " not defined because of singularities)\n", sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if(!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4, dimnames=list(cn, colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits=digits, signif.stars=signif.stars, na.print="NA", ...)
    }
    ##
    cat("\nResidual standard error:",
	format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    if (!is.null(x$fstatistic)) {
        fpval <- if(length(x$fstatistic)==3)
            format.pval(pf(x$fstatistic[1L], x$fstatistic[2L],
                           x$fstatistic[3L], lower.tail = FALSE), digits=digits)
        else
            format.pval(pchisq(x$fstatistic[1L], x$fstatistic[2L],
                               lower.tail = FALSE), digits=digits)
        
	cat("Multiple R-squared:", formatC(x$r.squared, digits=digits))
	cat(",\tAdjusted R-squared:",formatC(x$adj.r.squared,digits=digits),
	    "\nF-statistic:", formatC(x$fstatistic[1], digits=digits),
	    "on", x$fstatistic[2], "and",
	    x$fstatistic[3], "DF,  p-value:",
            fpval,
	    "\n")
    }
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



stata.regression <- function(x,
                             type = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"),
                             full = TRUE, mask,...) {
    type <- match.arg(type)
    
    
    sumx <- summary(x, type = type)
    Cf   <- coef(sumx)
    ci   <- confint(x, type = type, ...)
    nCf  <- nrow(Cf)
    vn   <- rownames(Cf)
    
    if(missing(mask))
        mask <- 1L:nCf

    Cf  <- Cf[mask,]
    nCf <- nrow(Cf)
    vn  <- vn[mask]
    ci  <- ci[mask,]
    vn <- abbreviate(vn, 12)
    
    mm <- match('(Intercept)', rownames(Cf))
    
    intercept <- 0
    if(!is.na(mm)) {
        rownames(Cf)[mm] <- '_cons'
        init <- 1
    }

    ddc <- cbind(ci, Cf)

    y <- names(x$model)[1]
    hh <- abbreviate(y, min = 11)

    if(nchar(hh)==11)
        substr(hh, start = 10, stop = 10) <- '~'   
    y <- format(hh, width = 12, justify = 'right')
    
    ina <- is.na(Cf[,1])
        
    if(intercept)
        ddc <- ddc[c(2:(nCf-1),1), c(3:5, 1:2)]
    else
        ddc <- ddc[, c(3:5, 1:2)]

    prettyfy <- function(x, width, digits, strip0 = FALSE, fix = FALSE) {
        if(fix)
            digit=digits
        else
            digit <- max(0, min(digits, digits-log10(max(abs(x)))))
        
        x <- formatC(x, format = 'f', width = width,
                     digits = digit,
                     preserve.width = 'common', flag="#")
        if(strip0) {
            x <- gsub("-0.", " -.", x, fixed = TRUE)
            x <- gsub("0.",   " .", x, fixed = TRUE)
        }
        x
    }

    if(type=='const')
        pv <- pt(abs(ddc[,3]), lower = FALSE, df = x$rank-1)*2
    else
        pv <- pnorm(abs(ddc[,3]), lower = FALSE)*2

    pv <- zapsmall(pv, digits = 4)
    fpv <- formatC(pv, 'f', width = 7, digits = 4, preserve.width = 'common', flag="#")
    oddc <- matrix(0, nrow(ddc), ncol(ddc))
    for(j in 1:nCf) {
        if(!is.na(ddc[j,1])) {
            oddc[j,1] <- prettyfy(ddc[j,1], 10, 7, TRUE)
            oddc[j,2] <- prettyfy(ddc[j,2], 10, 7, TRUE)
            oddc[j,3] <- prettyfy(ddc[j,3], 8, 2,  FALSE, TRUE)
            oddc[j,4] <- prettyfy(ddc[j,4], 10, 7, TRUE)
            oddc[j,5] <- prettyfy(ddc[j,5], 10, 7, TRUE)
        }
    }
            
    ddc <- oddc


    ## fCf <- formatC(ddc[,1:2], 'f', width = 10, digits = max(0, 7-log10(max(abs(ddc[!ina,1:2])))),
    ##                preserve.width = 'common', flag="#")
    ## fCf <- apply(fCf, 2, function(u) gsub("-0.", " -.", u, fixed = TRUE))
    ## fCf <- apply(fCf, 2, function(u) gsub("0.", " .", u, fixed = TRUE))

    ## fts <- formatC(ddc[,3, drop = FALSE], 'f',  width = 8, digits = max(0, 3-log10(max(abs(ddc[!ina,3])))),
    ##                preserve.width = 'common', flag="#")
    
    ## fCi <- formatC(ddc[,4:5], 'f', width = 10, digits = max(0, 7-log10(max(abs(ddc[!ina,4:5])))),
    ##                preserve.width = 'common', flag="#")
    ## fCi <- apply(fCi, 2, function(u) gsub("-0.", " -.", u, fixed = TRUE))
    ## fCi <- apply(fCi, 2, function(u) gsub("0.", " .", u, fixed = TRUE))


    xn <- cbind(ddc[,1:3], fpv, ddc[,4:5])
    vn <- formatC(vn, width = 12, flag = "#")


#    xn <- rownames(ddc)

#    Cf <- printMatCoef(ddc,  cs.ind = c(1:4)) 
#    Cf[Cf == '< 2.2e-16'] <- '0.0000'
#    Cf[which(as.numeric(Cf[,6])<1e-04),6] = '0.0000'
#    coeff <- format(ddc[,1], width = 10, digits = 5)
    
    if(full) {
        vv  <- matrix(0,1,6)
        vv[1,1] <- sumx$r.squared
        vv[1,2] <- sumx$adj.r.squared
        vv[1,3] <- sqrt(sum(residuals(sumx)^2)/x$df)
        vv[1,4] <- if(length(x$fstatistic)==3)
            pf(sumx$fstatistic[1L], sumx$fstatistic[2L], sumx$fstatistic[3L], lower.tail = FALSE)
        else
            pchisq(sumx$fstatistic[1L], sumx$fstatistic[2L], lower.tail = FALSE)
        
        rank <- formatC((x$rank-1), 'd', width = 3, digits = 0, flag = "# ")
        dfr  <- formatC(x$df+x$rank, 'd', width = 6, digits = 0, flag = "# ")
        vv <- zapsmall(vv, 5)

        fvv <- formatC(vv, format = 'f',
                       width = 8, digits = 4,
                       flag = "#")

        Fs <- formatC(sumx$fstatistic[1], 'f', width = 9, digits = max(0, 4-log10(max(abs(sumx$fstatistic[1])))), flag = "#")
        
        cat('Linear regression                                    Number of obs = ', dfr,'\n')
        cat('                                                     F(',rank,',', dfr,') = ', Fs,'\n', sep='')
        cat('                                                     Prob > F      = ', fvv[1,4],'\n')
        cat('                                                     R-squared     = ', fvv[1,1],'\n')
        cat('                                                     Adj R-squared = ', fvv[1,2],'\n')
        cat('                                                     Root MSE      = ', fvv[1,3],'\n')
    }

    colnames(xn) <- rep('',6)
    cat('------------------------------------------------------------------------------\n')
    if(type!="const")
        cat('                            Robust\n')
    cat(y,' |     Coef.   Std. Err.      t     P>|t|     [95% Conf. Interval]\n', sep = '')
    cat('------------------------------------------------------------------------------\n')
    for(j in 1:nCf)
        cat(vn[j],' |', xn[j,1],' ',xn[j,2],' ', xn[j,3],'  ',xn[j,4],'   ',xn[j,5],'  ',xn[j,6],'\n', sep ='')

    cat('------------------------------------------------------------------------------\n')
}


stripzero <- function(string)
{
    if(is.matrix(string))
        string <- apply(str, 1, function(x) if(!is.na(pmatch('0.',x))) paste('.',strsplit(x, '0.')[[1]][2], sep =''))
    if(is.array(string)| length(string)>1)
        for(j in 1:length(string))
            if(!is.na(pmatch('0.',string[j])))
                string[j] <- paste('.',strsplit(string[j], '0.')[[1]][2], sep = '')
    if(!is.na(pmatch('0.',string)))
        string <- paste('.',strsplit(string, '0.')[[1]][2], sep = '')
    return(string)
}
    

eviews.regression <-
    function(x, robust=TRUE, out='latex',...) 
{
    require(lmtest, quietly = TRUE)
    attx <- attributes(x)
    ff <- attx$factors
    yv <- rownames(ff)[ff==0]
    cat("\nDependent Variable: ", yv, "\nMethod: Least Squares\n")
    cat("Date:", format(Sys.time(), "%m-%d-%Y"),
        '  Time:', format(Sys.time(), "%H:%M:%S"),'\n')
    cat("Included Observations: ", length(x$fitted.values), ' after adjustments\n')
    if(!is.null(attx$iid) && !attx$iid)
        cat('White Heteroskedasticity-Consistent Standard Errors & Covariance\n')
    cat('======================================================================\n')
    cat('                    Coefficients   Std.Error   t-Statistic     Prob.\n')
    cat('======================================================================\n')
        
    cout <- coef(x)
    if(robust)
        require(sandwich)
    sx <- summary(x, robust=robust )
    resid <- x$resid
    ssq <- sum(resid^2)
    n <- length(resid)
    k <- length(coef(x))
    mX <- model.matrix(x)


    OUTbelow <- matrix(0,14,1)

    OUTabove <- coef(sx)
    if(rownames(OUTabove)[1]=='(Intercept)')
        rownames(OUTabove)[1]='C                 '
    colnames(OUTabove) <- rep("",4)
    colnames(OUTbelow) <- ""
    
    OUTbelow[1,] <- sx$r.squared
    OUTbelow[2,] <- sx$adj.r.squared
    OUTbelow[3,] <- sqrt(ssq/(n-k))
    OUTbelow[4,] <- ssq
    OUTbelow[5,] <- logLik(x)
    OUTbelow[6,] <- sx$fstatistic[1]
    OUTbelow[7,] <- pf(sx$fstatistic[[1]], df1=sx$fstatistic[[2]],
                       df2=sx$fstatistic[[3]], lower=FALSE)
    OUTbelow[8,] <- mean(x$effects)
    OUTbelow[9,] <- sd(x$effects)
    OUTbelow[10,] <- AIC(x)
    OUTbelow[11,] <- AIC(x, k=log(n)) ## SBIC
    OUTbelow[12,] <- AIC(x, k=((k/2)*log(n) + (1/2)*log(det(crossprod(mX))))/k)
    OUTbelow[13,] <- dwtest(x)$stat
    OUTbelow[14,] <-  NaN
    print(formatC(OUTabove,digits=6, format = 'f', width = 12), quote = FALSE, digits = 6)
    cat('======================================================================\n')
    no <- c('R-squared         ',
            'Adjusted R-squared',
            'S.E. of regression',
            'Sum squared resid ',
            'Log likelihood    ',
            'F-statistic       ',
            'Prob(F-statistic) ',
            '  Mean dependent var     ',
            '  S.D. dependent var     ',
            '  Akaike info criterion  ',
            '  Schwarz criterion      ',
            '  Hannan-Quinn criter.   ',
            '  Durbin-Watson stat     ',
            '                         ')

    for (j in 1:6)
        cat(no[j],
            formatC(OUTbelow[j,1], digits = 6, format = 'f',
                    width = 12),
            no[j+7],
            formatC(OUTbelow[j+7,1], digits = 6, format = 'f',
                    width = 12), '\n')
            
        cat(no[7],
            formatC(OUTbelow[7,1], digits = 6, format = 'f',
                    width = 12), '\n')
            cat('======================================================================\n\n')

    
}

lmlatexline <- function(x, ...)
    UseMethod('lmlatexline')

lmlatexline.reg <- function(object, vcov., ..., inline = FALSE, 
                            digits = max(3, getOption("digits") - 3),
                            scipen = options('scipen')[[1]], dmath = FALSE, r2 = FALSE) {
    coeff <- coef(object)
    if(missing(vcov.))
        vcov. <- vcov.reg
    se <- sqrt(diag(vcov.(object, ...)))
    sc <- ifelse(coeff>0,"+","-")
    cf <- names(coeff)
    ff <- attributes(object$terms)$factors
    yv <- rownames(ff)[rowSums(ff) == 0]
    trio <- format(round(cbind(abs(coeff), se), digits), digits = digits, trim = TRUE,
                   scientific = scipen)
    if(r2)
        r2 <- summary(object)$r.squared
    else
        r2 <- FALSE
    lmlatexline.print(trio, sc, yv, cf, dmath, r2, inline)
}


lmlatexline.ivreg <- function(object, vcov., ..., inline = FALSE, 
                            digits = max(3, getOption("digits") - 3),
                            scipen = options('scipen')[[1]], dmath = FALSE) {
    coeff <- coef(object)
    if(missing(vcov.))
        vcov. <- vcov.ivreg
    se <- vcov.(object, ...)
    sc <- ifelse(coeff>0,"+","-")
    cf <- names(coeff)
    ff <- attributes(object$terms)$factors
    yv <- rownames(ff)[rowSums(ff) == 0]
    trio <- format(round(cbind(abs(coeff), se), digits), digits = digits, trim = TRUE,
                   scientific = scipen)
    r2 <- FALSE
    lmlatexline.print(trio, sc, yv, cf, dmath, r2, inline)
}

lmlatexline.print <- function(trio, sc, yv, cf, dmath, r2, inline) {
    pstd <- function(z)
        paste('\\underset{(',z[2],')}{',z[1],'}\\times ',rownames(z), sep = '')

    if(cf[[1]] == '(Intercept)')
        rs <- paste(yv, "=", "\\underset{(", trio[1,2], ")}{", trio[1,1], "}", sep = '')
    else
        rs <- paste(yv,'=', sc[1], pstd(trio[1,,false]))
        
    for (j in 2:NROW(trio))
        rs <- paste(rs, sc[j], pstd(trio[j,,drop = FALSE]), paste = '')

    if(r2)
    {
        if(dmath)
            rs <- paste(rs, ',\\condition[]{$R^2 =', format(r2, digits = 3),'$}', sep = '')
        else
            rs <- paste(rs, ',\\quad R^2 =', format(r2, digits = 3), sep = '')
    }

    if(inline)
        cat("\n$",rs,"$\n")
    else
    {
        if(dmath)
            cat("\n\\begin{dmath*}")
        else
            cat("\n\\begin{equation*}")
        cat("\n",rs,"\n")
        if(dmath)
            cat("\\end{dmath*}\n\n")
        else
            cat("\\end{equation*}\n\n")
    }
}

vcov.reg <- function (object,
                      type = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"),
                      ...) {

    iscluster <- !is.null(object$cluster)
    type <- match.arg(type)
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
    V      <- matrix(NA, length(Qr$pivot), length(Qr$pivot))
    if(!iscluster) {
        if(type!="HAC") {
            meat.  <- meatHC(object, type = type, ...)            
            meat.  <- meat.[nNA, nNA]
            bread. <- bread(object, type = type, ...)
            Vout   <- (1/nrow(estfun(object)) * (bread. %*% meat. %*% bread.))            
        }
        else {
            Vout <- vcovHAC(object,  ...)
        }
    } else {
        r    <- z$residuals
        f    <- z$fitted.values
        w    <- if(is.null(z$weights)) rep(1,n) else z$weights
        facj <- n/(n-p)
        
        X          <- model.matrix(z)
        
        cluster    <- as.factor(z$cluster)
        nc         <- length(levels(cluster))
        j          <- order(cluster)
        clus.size  <- table(cluster)
        clus.start <- c(1, 1 + cumsum(clus.size))
        ## The dimension of X needs to be adjusted 
        X <- X[j, nNA, drop = FALSE]        
        w <- w[j]
        r <- r[j]    
        r <- r*w
        ##R <- solve(crossprod(X*sqrt(w)))
        R <- chol2inv(Qr$qr[p1, p1, drop = FALSE]) ## solve(X'X)
        
        mr2 <- function() {
            res <- NULL
            for (jj in 1:nc) {
                ind <- clus.start[jj]+ (0:(clus.size[jj]-1)) 
                Xi  <- X[ind,,drop=FALSE]
                Hgg <- chol(diag(length(ind))-Xi%*%R%*%t(Xi), pivot = TRUE)
                pivot <- attr(Hgg, "pivot")
                oo <- order(pivot)
                Hgg <- Hgg[,oo]
                res <- c(res, solve(Hgg)%*%r[ind])
            }
            res
        }
        
        mr3 <- function() {
            res <- NULL
            for (jj in 1:nc) {
                ind <- ind <- clus.start[jj]+ (0:(clus.size[jj]-1)) 
                Xi  <- X[ind,,drop=FALSE]
                Hgg <- solve(diag(length(ind))-Xi%*%R%*% t(Xi))
                res <- c(res, Hgg%*%r[ind])
            }
            sqrt((nc-1)/nc)*res
        }
        
        res <- switch(EXPR = type,                       
                      HC2 = mr2(),
                      HC3 = mr3(),
                      sqrt((n-1)/(n-p) * nc/(nc-1))*r
                      )
        
        score <- X*c(res)
        clus.start <- clus.start[-(nc + 1)]
        storage.mode(clus.start) <- "integer"
                                        #storage.mode(score) <- "double"
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

print.reg.wb <- function(x, ...)
{
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

summary.reg.wb <- function(x, ...)
{
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



print.summary.reg.wb <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    b <- coef(x)
    nb <- names(b)
    k <- length(b)
##  print(summary.reg(x))
##  NOTING AT THE MOMENT
    return()
##    cat('\nNull Hypothesis:')
##     for(j in 1:k)
##     {
##         cat('\n ', nb[j], ' = ', x$null[j], '\n')
##         cat(' wb p-value (a):' , format.pval(x$pv.tstat0[j]),
##             symnum(x$pv.tstat0[j], corr = FALSE, na = FALSE, 
##                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
##                    symbols = c("***", "**", "*", ".", " ")), '\n')
##         cat(' wb p-value (b):' , format.pval(x$pv.tstat0.se[j]),
##             symnum(x$pv.tstat0.se[j], corr = FALSE, na = FALSE, 
##                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
##                    symbols = c("***", "**", "*", ".", " ")), '\n')
##         cat(' wb p-value (c):' , format.pval(x$pv.tstat[j]),
##             symnum(x$pv.tstat[j], corr = FALSE, na = FALSE, 
##                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
##                    symbols = c("***", "**", "*", ".", " ")), '\n')
##         cat(' wb p-value (d):' , format.pval(x$pv.tstat.se[j]),
##             symnum(x$pv.tstat.se[j], corr = FALSE, na = FALSE, 
##                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
##                    symbols = c("***", "**", "*", ".", " ")), '\n')
##     }


##     cat("\n ---\n Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
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

format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    ## Manually exported!!
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
	  "%")

confint.reg <- function (object, parm, level = 0.95,
                         type = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"), ...) 
{
    ## ... other argument passsed down to vcov.reg (for instance when type is HAC
    type <- match.arg(type)
    cf <- coef(object)
    pnames <- names(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- format.perc(a, 3)
    fac <- qnorm(a)
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
    
fswv <- function(x, digits = 3)
    formatC(x, format = 'f', digits = digits)




printMatCoef <-
    function(x, digits = max(3, getOption("digits") - 2),
	     signif.stars = getOption("show.signif.stars"),
             signif.legend = signif.stars,
	     dig.tst = max(1, min(5, digits - 1)),
	     cs.ind = 1:k, tst.ind = k+1, zap.ind = integer(0),
	     P.values = NULL,
	     has.Pvalue = nc >= 4 && substr(colnames(x)[nc], 1, 3) == "Pr(",
             eps.Pvalue = .Machine$double.eps,
	     na.print = "NA", ...)
{
    ## For printing ``coefficient matrices'' as they are in summary.xxx(.) where
    ## xxx in {lm, glm, aov, ..}. (Note: summary.aov(.) gives a class "anova").

    ## By Default
    ## Assume: x is a matrix-like numeric object.
    ## ------  with *last* column = P-values  --iff-- P.values (== TRUE)
    ##	  columns {cs.ind}= numbers, such as coefficients & std.err  [def.: 1L:k]
    ##	  columns {tst.ind}= test-statistics (as "z", "t", or "F")  [def.: k+1]

    if(is.null(d <- dim(x)) || length(d) != 2)
	stop("'x' must be coefficient matrix/data frame")
    nc <- d[2L]
    if(is.null(P.values)) {
        scp <- getOption("show.coef.Pvalues")
        if(!is.logical(scp) || is.na(scp)) {
            warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
            scp <- TRUE
        }
	P.values <- has.Pvalue && scp
    }
    else if(P.values && !has.Pvalue)
	stop("'P.values' is TRUE, but 'has.Pvalue' is not")

    if(has.Pvalue && !P.values) {# P values are there, but not wanted
	d <- dim(xm <- data.matrix(x[,-nc , drop = FALSE]))
	nc <- nc - 1
	has.Pvalue <- FALSE
    } else xm <- data.matrix(x)

    k <- nc - has.Pvalue - (if(missing(tst.ind)) 1 else length(tst.ind))
    if(!missing(cs.ind) && length(cs.ind) > k) stop("wrong k / cs.ind")

    Cf <- array("", dim=d, dimnames = dimnames(xm))
browser()
    ok <- !(ina <- is.na(xm))
    ## zap before deciding any formats
    for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
    if(length(cs.ind)) {
	acs <- abs(coef.se <- xm[, cs.ind, drop=FALSE])# = abs(coef. , stderr)
	if(any(ia <- is.finite(acs))) {
	    ## #{digits} BEFORE decimal point -- for min/max. value:
	    digmin <- 1 + if(length(acs <- acs[ia & acs != 0]))
		floor(log10(range(acs[acs != 0], finite = TRUE))) else 0
            Cf[,cs.ind] <- format(round(coef.se, max(1,digits-digmin)),
                                  digits=digits)
        }
    }
    
    if(length(tst.ind))
	Cf[, tst.ind]<- format(round(xm[, tst.ind], digits = dig.tst),
                               digits = digits)
    if(any(r.ind <- !((1L:nc) %in%
                      c(cs.ind, tst.ind, if(has.Pvalue) nc))))
	for(i in which(r.ind)) Cf[, i] <- format(xm[, i], digits=digits)
    okP <- if(has.Pvalue) ok[, -nc] else ok
    ## we need to find out where Cf is zero.  We can't use as.numeric
    ## directly as OutDec could have been set.
    ## x0 <- (xm[okP]==0) != (as.numeric(Cf[okP])==0)
    x1 <- Cf[okP]
    dec <- getOption("OutDec")
    if(dec != ".") x1 <- chartr(dec, ".", x1)
    x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
    if(length(not.both.0 <- which(x0 & !is.na(x0)))) {
	## not.both.0==TRUE:  xm !=0, but Cf[] is: --> fix these:
	Cf[okP][not.both.0] <- format(xm[okP][not.both.0],
                                      digits = max(1, digits-1))
    }
    if(any(ina)) Cf[ina] <- na.print
    if(P.values) {
        if(!is.logical(signif.stars) || is.na(signif.stars)) {
            warning("option \"show.signif.stars\" is invalid: assuming TRUE")
            signif.stars <- TRUE
        }
	if(any(okP <- ok[,nc])) {
	pv <- as.vector(xm[, nc]) # drop names
	    Cf[okP, nc] <- format.pval(pv[okP],
                                       digits = dig.tst, eps = eps.Pvalue)
	    signif.stars <- signif.stars && any(pv[okP] < .1)
	    if(signif.stars) {
		Signif <- symnum(pv, corr = FALSE, na = FALSE,
				 cutpoints = c(0,  .001,.01,.05, .1, 1),
				 symbols   =  c("***","**","*","."," "))
		Cf <- cbind(Cf, format(Signif)) #format.ch: right=TRUE
	    }
	} else signif.stars <- FALSE
    } else signif.stars <- FALSE
    #print.default(Cf, quote = FALSE, right = TRUE, na.print=na.print, ...)
    #if(signif.stars && signif.legend)
    #    cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")
    invisible(Cf)
}

