##' @export
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

    xn <- cbind(ddc[,1:3], fpv, ddc[,4:5])
    vn <- formatC(vn, width = 12, flag = "#")
    
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

##' @export
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

##' Output regression in latex form on a line
##'
##' Details
##' @title lmlatexline
##' @param x  object
##' @param ... other arguments
##' @rdname lmlatexline
##' @export 
lmlatexline <- function(x, ...)
    UseMethod('lmlatexline')

##' @S3method lmlatexline reg
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

##' @S3method lmlatexline ivreg
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

format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    ## Manually exported!!
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
	  "%")

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
				 cutpoints = c(0, .001,.01,.05, .1, 1),
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
