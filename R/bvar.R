##################################################################
## bvar.R /
## Author: Giuseppe  Ragusa 
## Time-stamp: "2010-04-07 19:17:27 gragusa bvar.R" 
##
## Description: Bayesian BVAR
##################################################################

lags <- function(x, k = 1)
{
    T <- NROW(x)
    nc <- NCOL(x)
    if(nc>1) {
        for(j in 1:nc)
            x[,j] <- c(rep(NA, k), x[-((T-k+1):T),j])
    }
    else
        x <- c(rep(NA, k), x[-((T-k+1):T)])
    x
}

vars <- function(Y, p = 1, h = 1, const = TRUE)
{
    ## p var lag
    ## h step head for direct forecast
    ## ==========================================
    ## Yt = Y(t-h) Y(t-h-1) Y(t-h-2) ... Y(t-h-p)
    ## ==========================================
    ## const = TRUE, intercept
    ##       = FALSE, no intercept
    ## 
    ## if(class(Y)[1]=='zooreg'|class(Y)[1]=='zoo')
    ##     time <- time(Y)

    ## Y is either a mts or zoo

    #if(class(Y)[1]!='mts' & class(Y)[1]!='zooreg' & class(Y)[1]!='zoo')
    #    stop()
    time <- index(Y)
    
    if(!is.numeric(p))
        stop("'p' must be an integer greater than 0")
    if(length(p)>1)
        stop("'p' must be an integer greater than 0")

    if(!is.numeric(h))
        stop("'h' must be an integer greater than 0")
    if(length(h)>1)
        stop("'h' must be an integer greater than 0")

    nc <- NCOL(Y)
    T  <- NROW(Y)
    ## Create Matrix of lags
    XX <- lags(Y, k = h)
    if(p>1)
        for(j in 2:p) 
            X <- XX <- cbind(XX, lags(Y, k = (h+j-1)))
    else
        X <- XX
    YY <- matrix(Y, ncol = NCOL(Y), nrow = NROW(Y))
    XX <- matrix(XX,ncol = NCOL(XX), nrow = NROW(XX))
    
    YYXX <- na.omit(YX <- cbind(YY, XX))
    YY <- (YYXX[,1:nc])
    XX <- (YYXX[,-(1:nc)])
    ncx <- NCOL(XX)
    nrx <- NROW(XX)
    if(const)
        XX <- cbind(1,XX)
    ## OLS Coefficients
    betas <- solve(cXX <- crossprod(XX), crossprod(XX,YY))

    vvars <- array(list(matrix(0,ncx,ncx)), nc) 
    ## OLS Variances
    for(j in 1:nc)
    {
        eta <- YY[,j]-tcrossprod(betas[,1], XX)
        sigmahat <- sum(eta^2)/(nrx-ncx)
        vvars[[j]] <- sigmahat*solve(cXX)
    }

    ##    Y <- as.zoo(YX[,1:nc])
    ##    X <- as.zoo(YX[,-(1:nc)])
    ##    time(Y) <- time
    ##    time(X) <- time
    
    out <- list(coefficients = betas, vcov = vvars,
                p = p, h = h,
                YY = YY, ## Matrix
                XX = XX, ## Matrix
                Y  = Y,   ## Original Data
                X  = X,   ## Original Data
                const = const)
    class(out) <- 'VAR'
    out
}


bvar <- function(Y, p = 1, h = 1,  
                 lambda0, lambda1, lambda2,
                 si, si.type = c('null', 'ar'),
                 c0 = 0,
                 b0.ownlag = rep(1, ncol(Y)),
                 b0.offdia = 0,
                 k.type = c('squared', 'linear','exponential'))
{
    ## lambda0 - constant
    ## lambda1 - own lag
    ## lambda2 - other lags
    si.type <- match.arg(si.type)
    k.type  <- match.arg(k.type)
    T <- NROW(Y)
    m <- NCOL(Y)

    if(length(c0)>1)
        stop()
    if(length(b0.ownlag)>m)
        stop()
    if(length(b0.offdia)>1)
        stop()
  
    constant <- ifelse(is.null(c0), FALSE, TRUE)
    cadj <- ifelse(is.null(c0), 0, 1)
    pm <- p*m+cadj
    out <- vars(Y, p = p, h = h, const = constant)
    ## YY is Txp (class mts) or (zoo)...
    ## XX is Txpm
    ##YY <- as.matrix(YY)
    ##XX <- as.matrix(out$X)
    ##if(constant)
    ##    XX <- cbind(1,XX)

    ##YYXX <- na.omit(cbind(YY,XX))
    ##YY <- YYXX[,1:m]
    ##XX <- YYXX[,-(1:m)]
    YY <- out$YY
    XX <- out$XX

    ## Construct B0
    B0 <- matrix(b0.offdia, p*m, m)
    diag(B0) <- b0.ownlag
    if(constant)
        B0 <- rbind(c0, B0)

    if(missing(si)) { ## Calculate si
        si <- array(0, m)
        if(si.type=='null') {
            u  <- YY-XX%*%B0  ## Tx(m+cadj)
            si <- apply(u, 2, sd)
        }
        else {
            for(j in 1:m) {
                fm <- arima(YY[,j], order = c(p, 0, 0))
                si[j] <- sqrt(fm$sigma2)
            }
        }}
    
    ## Prior Variances
    ## X
    ##
    ## Xt-1 ... Xt-1 Xt-2 .. Xt-2
    ## Own lags are at:
    ## (1,m+1, 2*m+1,...)
    
    Sigma0 <- matrix(0, p*m, m)
    
    for(i in 1:m) {
        ol <- seq(i,p*m,by=m)   ## Own lags
        jl <- 1:(p*m)
        jl <- jl[-match(ol,jl)]
        kol <- 1:length(ol)     ## Own lag k
        kjl <- rep(1:p, each=(m-1))  ## Other lags
        if(k.type=='squared') {
            kol <- kol^2
            kjl <- kjl^2
        }
        if(k.type=='exponential') {
            kol <- exp(kol)
            kjl <- exp(kjl)
        }
            
        idx <- rep(jl[1:(m-1)], p)
        
        Sigma0[ol,i] <- lambda1/kol
        ls <- lambda2*si[i]^2
        ls <- ls/si[idx]^2
        Sigma0[jl,i] <- ls/kjl
    }

    if(constant) {
        S00 <- lambda0*si^2
        Sigma0 <- rbind(S00, Sigma0)
    }
    
    ## Calculate Posterior Moments

    gammai <- matrix(0, pm, m)
    Sigmab <-  array(list(matrix(0, (m+cadj),(m+cadj))), m) 
    for(j in 1:m)
    {
        Sigmat <- diag(Sigma0[,j]^(-1))
        Sigmab[[j]] <- 
            solve(Sigmat+si[j]^(-2)*crossprod(XX))
        gammai[,j] <-
            Sigmab[[j]]%*%(Sigmat%*%B0[,j]+si[j]^(-2)*crossprod(XX,YY[,j]))
    }
    out <- list(coefficients = gammai, vcov = Sigmab,
                p = p, h = h,
                YY        = out$YY,
                XX        = out$XX,
                Y         = out$Y,
                X         = out$X,
                lambda0   = lambda0,
                lambda1   = lambda1,
                lambda2   = lambda2,
                si        = si,
                si.type   = si.type,
                c0        = c0,
                b0.ownlag = b0.ownlag,
                b0.offdia = b0.offdia,
                k.type    = k.type,
                const     = constant,
                OLS       = out)
    class(out) <- c('VAR', 'BVAR')
    out
}


rw.rolling.forecast <- function(Y, h)
{
    nvar <- NCOL(Y)
    T    <- NROW(Y)

    forecast       <- matrix(NA, T, nvar)
    forecast       <- as.zoo(forecast)
    time(forecast) <- time(Y)

    forecast <- lags(Y, k = n.ahead)
    out = list(forecast = forecast, to.forecast = Y, n.ahead = n.ahead)
    class(out) <- c('rolling.forecast')
    out
}

predict.out.vars <- function(obj, newdata)
{
    
    ## New data is passed without intercept
    if(obj$const)
        crossprod(coef(obj), c(1,newdata))
    else
        crossprod(coef(coef), newdata)
}

bvar.rolling.forecast <- function(obj, window.size = 50, h)
{
    rw <- window.size
    p <- obj$p
    if(missing(h))
        h <- obj$h
    Y <- obj$Y
    X <- obj$X
    nvar <- NCOL(Y)
    T <- NROW(Y)
    
    forecast <- matrix(NA, T, nvar)
    forecast <- as.zoo(forecast)
    time(forecast) <- time(Y)
    fq <- frequency(forecast)

    Y <- as.zoo(Y)
    X <- as.zoo(obj$X)
    colnames(X) <- NULL
    dtx <- time(X)
    ##X <- as.matrix(X)
    for(j in 1:(T-rw-(h-1)))
    {
        ## Estimate
        ## Yt = Y(t-h-1) Y(t-h-2) ... Y(t-h-p)
        Ytmp <- Y[(1+j-1):(rw+j-1), ,drop = FALSE]
        vr  <- bvar(Ytmp, p   = p, h = h,
                    lambda0   = obj$lambda0,
                    lambda1   = obj$lambda1,
                    lambda2   = obj$lambda2,
                    si        = obj$si,
                    si.type   = obj$si.type,
                    c0        = obj$c0,
                    b0.ownlag = obj$b0.ownlag,
                    b0.offdia = obj$b0.offdia,
                    k.type    = obj$k.type)
        dty <- time(Ytmp[rw,])+1/fq;
        mtc <- match(dty, dtx)
        nvl <- X[mtc,]
        forecast[(rw+j+h-1),] <- predict.out.vars(vr, nvl)
    }
    fc  <- as.zoo(forecast)
    fq <- frequency(forecast)
    time(forecast) <- time(forecast)
    out <- list(forecast = fc, to.forecast = Y, h = h)
    class(out) <- c('rolling.forecast')
    out
}

########
var.rolling.forecast <- function(obj, window.size = 50, h)
{
    rw <- window.size
    p <- obj$p
    if(missing(h))
        h <- obj$h
    Y <- obj$Y
    X <- obj$X
    nvar <- NCOL(Y)
    T <- NROW(Y)
    
    forecast <- matrix(NA, T, nvar)
    forecast <- as.zoo(forecast)
    time(forecast) <- time(Y)
    fq <- frequency(forecast)

    Y <- as.zoo(Y)
    X <- as.zoo(obj$X)
    colnames(X) <- NULL
    dtx <- time(X)
    ##X <- as.matrix(X)
    for(j in 1:(T-rw-(h-1)))
    {
        ## Estimate
        ## Yt = Y(t-h-1) Y(t-h-2) ... Y(t-h-p)
        Ytmp <- Y[(1+j-1):(rw+j-1), ,drop = FALSE]
        vr  <- vars(Ytmp, p = p, h = h, const = obj$const)
        dty <- time(Ytmp[rw,])+1/fq;
        mtc <- match(dty, dtx)
        nvl <- X[mtc,]
        forecast[(rw+j+h-1),] <- predict.out.vars(vr, nvl)
    }
    fc  <- as.zoo(forecast)
    fq <- frequency(forecast)
    time(forecast) <- time(forecast)
    out <- list(forecast = fc, to.forecast = Y, h = h)
    class(out) <- c('rolling.forecast')
    out
}

########

summary.rolling.forecast <- function(obj)
{
    nc    <- NCOL(obj$forecast)
    cln   <- na.omit(cbind(obj$forecast, obj$to.forecast))
    cast  <- cln[,1:nc]
    vcast <- cln[,-(1:nc)]
    
    mse <-
        apply(vcast-cast, 2,
              function(x) sqrt(mean(sum(x^2))))
    mae <-
        apply(vcast-cast, 2,
              function(x) mean(abs(x)))

    summ <- summary(vcast-cast)
    #print(summ)
    #cat('\n')
    #cat('MSE: ', mse, 'MAE: ', mae,'\n')
    list(mse = mse, mae = mae, summary = summ)
}
              
plot.rolling.forecast <- function(x, which = 'all',...,z)
{
    nc    <- NCOL(x$forecast)
    vn    <- colnames(x$to.forecast)
    cln   <- na.omit(cbind(x$forecast, x$to.forecast))
    it    <- time(na.omit(x$forecast))
    cast  <- as.zoo(cln[,1:nc])
    vcast <- as.zoo(cln[,-(1:nc)])

    time(vcast) <- it
    time(cast) <- it

    if(is.character(which))
        if(which=='all')
            which <- 1:nc
    else
        if(!all(which<nc))
            stop("'which' must be a vector of indexes or 'all'")
    
    for(j in which)
    {
        aser <- c(as.numeric(vcast[,j]),as.numeric(cast[,j]))
        ylim <- c(min(aser), max(aser))
        plot(vcast[,j, drop = FALSE], ylab = '',
             xlab = 'T', type = 'l', main = vn[j], ylim = ylim, lwd = 2, ...)
        lines(cast[,j], type = 'l', col = 'gray', lwd = 2, lty = 1)
        if(!missing(z))
            lines(z[,j], type = 'l', col = 'darkgray', lwd = 2, lty = 1)
    }
}
