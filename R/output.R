##' "Statify" regression output
##'
##' Given a regression object out a stata-like table of results.
##' 
##' @title statify regression object
##' @param x a suitable object
##' @param vcov. the variance to be used
##' @param full if full output is needed
##' @param mask coefficient to be exluded
##' @param df degrees of freedom for the p-values
##' @param ... 
##' @return Output 
##' @author Giuseppe M. Ragusa
##' @export
statify <- function(x,
                    vcov. = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"),
                    full = TRUE, mask, df=NULL, ...) {
  if(missing(mask))
    mask <- 1L:length(coef(x))
  scipen <- options()$scipen
  type <- match.arg(vcov.) 
  coef <- coeftest(x, df = df)
  coef <- cbind(coef, confint(x, vcov. = type))  
  coef[,4] <- zapsmall(coef[,4], digits = 4)
  coef <- coef[mask,] 
  vn <- rownames(coef)
  mm <- match('(Intercept)', vn)
  if(!is.na(mm)) {
    vn[mm] <- '_cons'   
    vn   <- vn[c(2:length(vn),1)]
    coef <- coef[c(2:length(vn),1),]
  }
  fcoef <- regformat(format.df(coef, numeric=FALSE), len = c(7,7,3,3,7,7),
                     strip.zero = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE),
                     width = c(11, 11, 9, 8, 13, 12))
  colnames(fcoef) <- colnames(coef)
  rownames(fcoef) <- rownames(coef)
  cf <- fcoef
  
  y <- names(x$model)[1]
  
  vn <- sapply(vn, FUN=function(u) paste(substr('            ', 1, 12-nchar(u)), 
                                         u, sep=''))  
  if(full) {
    sumx <- summary(x)
    fstat <- sumx$fstatistic
    info <- matrix(0,5,1)
    info[1,1] <- fstat[1L]
    info[2,1] <- ifelse(length(fstat)==3, 
                        pf(fstat[1L], fstat[2L], fstat[3L], lower.tail = FALSE),
                        pchisq(fstat[1L], fstat[2L], lower.tail = FALSE))      
    info[2,1] <- zapsmall(info[2,1], digits=4)
    info[3,1] <- sumx$r.squared
    info[4,1] <- sumx$adj.r.squared
    info[5,1] <- sqrt(sum(residuals(sumx)^2)/x$df)
    info <- matrix(apply(info, 1, format, scientific=0), 1,5, byrow=TRUE)
    info <- regformat(info, len = c(4, 5, 5, 5, 6), strip.zero=rep(FALSE, 5),
                      width=rep(7,5))
    ind <- which(as.numeric(info)>10^5)
    if(length(ind))
      info[1,ind] <- paste(format.df(as.numeric(info[1,ind])/10^6, numeric=FALSE, dec=1), 'e+06', sep='')
    
    rank <- formatC((x$rank-1), 'd', width = 3, digits = 0, flag = "# ")
    dfr  <- format(ceiling(x$df+x$rank), scientific=(1-scipen))
    dfr  <- paste(strtrim('      ',6-nchar(dfr)), dfr, sep='')
    obs  <- format(ceiling(x$df), scientific=(1-scipen))
    obs  <- paste(strtrim('      ',6-nchar(dfr)), dfr, sep='')
    
    cat('Linear regression                                     Number of obs = ', dfr,'\n')
    cat('                                                      F(',rank,',', dfr,') = ', info[1,1],'\n', sep='')
    cat('                                                      Prob > F      =', info[1,2],'\n')
    cat('                                                      R-squared     =', info[1,3],'\n')
    cat('                                                      Adj R-squared =', info[1,4],'\n')
    cat('                                                      Root MSE      =', info[1,5],'\n')
  }
  
  cat('------------------------------------------------------------------------------\n')
  if(type!="const")
    cat('                             Robust\n')
  cat('           ',y,' |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]\n', sep = '')
  cat('------------------------------------------------------------------------------\n')
  for(j in 1:length(coef(x)))
    ##cat(vn[j],' |', cf[j,1],' ',cf[j,2],' ', cf[j,3],'  ',cf[j,4],'   ',cf[j,5],'  ',cf[j,6],'\n', sep ='')
    cat(paste(vn[j], ' |', cf[j,1], cf[j,2], cf[j,3], cf[j,4], cf[j,5], cf[j,6], sep=''),'\n')
  
  
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
lmlatexline.reg <- function(x, se= TRUE, vcov., ..., r2 = FALSE, 
                            dmath = FALSE, inline = FALSE, 
                            purge.factor.name = TRUE, purge.I = TRUE,
                            digits = max(3, getOption("digits") - 3),
                            scipen = options('scipen')[[1]]) 
  {
  coeff <- coef(x)
  if(missing(vcov.))
    vcov. <- vcov.reg
  serr <- sqrt(diag(vcov.(x, ...)))
  sc <- ifelse(coeff>0,"+","-")
  cf <- names(coeff)
  ff <- attributes(x$terms)$factors
  yv <- rownames(ff)[rowSums(ff) == 0]

  tmp <- cbind(abs(coeff), serr)
  ##tmp <- round(tmp_original, 2)
  tmp1 <- round(tmp[,1], digits = digits)
  tmp2 <- round(tmp[,2], digits = digits)
  
  digits0 <- digits
  while(any(tmp2==0)) {
    test <- tmp2==0
    digits0 <- digits0+1
    tmp2[which(test)] <- round(tmp[which(test),2], digits = digits0)      
  }
  
  digits0 <- digits
  while(any(tmp1==0)) {
    test <- tmp1==0
    digits0 <- digits0+1
    tmp1[which(test)] <- round(tmp[which(test),1], digits = digits0)      
  }
  
  trio <- format(cbind(tmp1,tmp2), digits=max(1, digits-1), drop0trailing=TRUE, scientific = 10000)
  if(nrow(trio)==1) rownames(trio) <- cf
  
  if(r2)
    r2 <- summary(x)$r.squared
  else
    r2 <- FALSE
  
  if(purge.factor.name) {
    topurge <- names(which(attributes(x$terms)$dataClasses=="factor"))
    for(j in topurge) {
      rownames(trio) <- str_replace_all(rownames(trio), j, "")   
    }
  }
  
  if(purge.I) {
    tmp <- str_locate(rownames(trio), "I\\(")
    whereI <- which(!is.na(tmp[,1]))
    for(j in whereI) {
      rownames(trio)[j] <- str_sub(rownames(trio)[j], start = 3, end = str_length(rownames(trio)[j])-1)
    }
  }
  
  
  
  lmlatexline.print(trio, sc, yv, cf, dmath, r2, se, inline)
}

##' @S3method lmlatexline lm
lmlatexline.lm <- lmlatexline.reg


##' @S3method lmlatexline ivreg
lmlatexline.ivreg <- function(x, se = TRUE., vcov., ..., 
                              dmath = FALSE, inline = FALSE, 
                              digits = max(3, getOption("digits") - 3),
                              scipen = options('scipen')[[1]]) {
  object <- x
  coeff <- coef(object)
  if(missing(vcov.))
    vcov. <- vcov.ivreg
  serr <- vcov.(object, ...)
  sc <- ifelse(coeff>0,"+","-")
  cf <- names(coeff)
  ff <- attributes(object$terms)$factors
  yv <- rownames(ff)[rowSums(ff) == 0]
  trio <- format(round(cbind(abs(coeff), serr), digits), digits = digits, trim = TRUE,
                 scientific = scipen)
  r2 <- FALSE
  lmlatexline.print(trio, sc, yv, cf, dmath, r2, se, inline)
}

lmlatexline.print <- function(trio, sc, yv, cf, dmath, r2, se, inline) {
  pstd <- function(z)
    paste('\\underset{(',z[2],')}{',z[1],'} \\,',rownames(z), sep = '')
  
  pstd_nose <- function(z) {
    paste(z[1], '\\,', rownames(z), sep = '')
  }
  
  if(cf[[1]] == '(Intercept)')
    if(se) {
      rs <- paste(yv, "=", "\\underset{(", trio[1,2], ")}{", trio[1,1], "}", sep = '')
    } else {
      rs <- paste(yv, "=",  trio[1,1], sep = '')
    }
  else
    if(se) {
      rs <- paste(yv,'=', sc[1], pstd(trio[1,,drop=FALSE]))
    } else {
      rs <- paste(yv,'=', sc[1], trio[1,1])
    }

  if(NROW(trio)>1) {
  for (j in 2:NROW(trio)) {
    if(se) {
      rs <- paste(rs, sc[j], pstd(trio[j,,drop = FALSE]), sep = '')
    } else {
      rs <- paste(rs, sc[j], pstd_nose(trio[j,,drop=FALSE]), sep="")
    }
  } 
}

  if(r2) {
    if(dmath)
      rs <- paste(rs, '\\condition[,]{$R^2 =', format(r2, digits = 3),'$}', sep = '')
    else
      rs <- paste(rs, ',\\quad R^2 =', format(r2, digits = 3), sep = '')
  }

  if(inline)
    cat("\n$",rs,"$\n")
  else {
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
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")

fswv <- function(x, digits = 3) {
  formatC(x, format = 'f', digits = digits)
}

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

    ## By Default:
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
	    Cf[okP, nc] <- base::format.pval(pv[okP],
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


stripzero <- function(string) {
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

regformat  <- function(x, len, strip.zero, width,...) {
  ## Stata Default for coef
  ##     len <- c(8,7,4,4,7,7)
  ##     strip.zero <- c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)
  nr  <- nrow(x)
  nc  <- ncol(x)  
  xout <- matrix('', ncol=nc, nrow=nr)
  for(h in 1:nc) {
    for(j in 1:nr){
      tmp <- strsplit(x[j,h], "\\.")[[1]]
      intg <- gsub("^\\s+|\\s+$", "", tmp[1])
      mant <- tmp[2]
      nintg <- nchar(intg)-ifelse(strip.zero[h], (intg=="0"|intg=="-0"), 0)
      nmant <- nchar(mant)
      decimal <- len[h]-nintg+(strsplit(intg, '-')[[1]][1]=='')
      if(decimal<=0){
        tmp <- format.df(ceiling(as.numeric(x[j,h])), numeric=FALSE)
      } else {
        tmp <- format.df(as.numeric(x[j,h]), dec = decimal, numeric.dollar=FALSE)
      }
      
      if(strip.zero[h]) {
        tmp0 <- strsplit(tmp, "\\.")[[1]]
        if(tmp0[1]=="0")
          tmp <- paste('.', tmp0[2], sep= '')
        if(tmp0[1]=="-0")
          tmp <- paste('-.', tmp0[2], sep= '')  
      }
      xout[j,h]  <- paste(substr('          ', 1, width[h]-nchar(tmp)), tmp, sep='')
    }  
  }
  ##matrix(format.df(as.numeric(xout), numeric=F), ncol=nc, nrow=nr)
  xout
}




