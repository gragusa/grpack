##################################################################
## jj.R /
## Author: Giuseppe  Ragusa 
## Time-stamp: "2010-07-04 20:17:47 gragusa mcmcMH.R" 
##
## Description: All purposes Metropolis Hasting algorithm
##################################################################


'mcmcMH' <- function(fun, theta.init, V, par.block,
                     mcmc = 10000, burnin = 10000, 
                     thin = 1, verbose = mcmc/10,
                     adaptive = FALSE, tune = 1,  rwmcmc = TRUE,
                     accept.target = 0.33, stop.adaptation = burnin,
                     gamma = 0.5, lbound = 10^-6, hbound = 100,
                     s0 = 1, propDist, debug = FALSE, ...)
{
    require('coda', quietly = TRUE)
    k  <- NROW(theta.init)
    nb <- 0
    if(!missing(par.block))
        if(is.list(par.block)) {
            nb <- length(par.block)
            ne <- length(upb <- unlist(par.block))
            if(ne!=length(theta.init))
                stop('parameters block is not compatible with theta.init')
            if(!all(sort(unlist(par.block))==c(1:k)))
                warning('not all indexes are updated') }
    
    if(missing(propDist))
        propDist <- if(rwmcmc)
            function(n) rnorm(n, mean = 0, sd = 1)
        else
            function(n) rt(n, df = 4)/sqrt(2)
    
    cholbV <- chol(V)
    total.simulations <- burnin + mcmc
    nsample <- round( mcmc/thin )
    
    count     <- 1
    ratioby   <- if (nb>0) array(0,nb) else 0
    accept    <- if (nb>0) array(0,k) else 0
    QQ        <- array(0, nsample, k)
    sample    <- matrix(0, nsample, k)
    theta_can <- theta <- theta.init
    
    Q <- fun(theta.init, ...) 
    if (Q == -sqrt(.Machine$double.xmax))
        stop('starting value must have >0 probability')
    
    QQ[1] <- Q
    for (iter in 2:(total.simulations+1)) {
        if (nb>0) {
            for (jj in 1:nb) {
                pb <- par.block[[jj]]
                theta_can <- theta
                if (rwmcmc)
                    theta_can[pb] <-
                        theta[pb] + propDist(length(pb))*tune[pb]*diag(V)[pb]^.5
                else
                    theta_can[pb] <- 
                        start+propDist(length(pb))*tune[pb]*diag(V)[pb]^.5
                
                Q_can <- fun(theta_can, ...)
                ratio <- exp(Q_can-Q)
                ratioby[pb] <- ratio
                if (runif(1) < ratio) {
                    theta <- theta_can
                    Q <- Q_can
                    accept[pb] <- accept[pb]+1
                }
                else
                    theta_can <- theta
            }
            QQ[iter] <- Q
        }
        else {
            if (rwmcmc) 
                theta_can <- 
                    theta+tune*t(cholbV)%*%propDist(k)
            else 
                theta_can <- 
                    theta.init+tune*t(cholbV)%*%propDist(k)
            
            Q_can <- fun(theta_can, ...) 
            ratio <- exp(Q_can-Q)
            
            if (runif(1) < ratio) {
                theta <- theta_can
                Q <- Q_can
                accept <- accept+1 }
        }
        if (iter%%thin == 0 && iter > burnin) {
            sample[ count, ] <- theta
            QQ[count] <- Q 
            count = count+1
        }
        if (is.numeric(adaptive) && adaptive > 0
            && iter%%adaptive == 0 && iter < stop.adaptation ) {
            if (nb>0)
                for (jjj in 1:nb) {
                    pb <- par.block[[jjj]]
                    tune[pb] <-
                        tune[pb] + (s0/iter^gamma)* (min(1,ratioby[jj])-accept.target)
                    tune[pb] <- ifelse(tune[pb] < lbound, lbound, tune[pb])
                    tune[pb] <- ifelse(tune[pb] > hbound, hbound, tune[pb])
                }
            else {
                cratio <- ifelse(is.finite(ratio),ratio,0)
                tune <- tune + (s0/iter^gamma)*( min(1,cratio) - accept.target )
                tune <- ifelse( tune < lbound, lbound, tune )
                tune <- ifelse( tune > hbound, hbound, tune ) }
        }
        
        if ( verbose > 0 && iter%%verbose == 0 ) {
            cat("---------------------------------------------\n")
            cat("Iteration: ", format(iter,width = 4)," Tune....:",format(tune,width = 4),"\n")
            if(nb>0)
                theta.out <-  cbind(theta, Q,  accept/iter)
            else
                theta.out <-  cbind(theta, rep(Q,k), rep(accept/iter))
            colnames(theta.out) <- c('Status','fun','Acceptance')
            print(theta.out,quote = FALSE,digits = 4,print.gap = 4,right = TRUE)
        }
    }
    
    names.theta <- names(start)
    
    if (is.null(names.theta) || length(names.theta) != k)
        names.theta <- paste("theta.",1:k,sep = "")
    
    colnames(sample)        <- c(names.theta)
    attr(sample,"tune")     <- tune
    attr(sample,"lik")      <- QQ
    attr(sample,'type')     <- ifelse(rwmcmc, 'random walk', 'independent')
    attr(sample,'proposal') <- propDist
    mcmc(data = sample, start = burnin + 1, thin = thin)
}



histogram.mcmc <- function(x, kernel.est.density = FALSE, normal.density = FALSE,
                           arg.mean, arg.sd, mask, xlab = "Posterior Distributions",
                           nint = round(log2(length(x)) + 1), col = "gray",
                           density.col = 'black', density.lty = 1, priors = NULL, ... )
{
    require(lattice)
    if (missing(mask))
        mask <- c(1:NCOL(x))
    colx <- NROW(mask)
    if (normal.density)
        summary.chain <- summary(x)
    if(missing(arg.mean) & normal.density)
        arg.mean <- matrix(summary.chain$statistics,NCOL(x),4)[mask,1]
    if(missing(arg.sd) & normal.density)
        arg.sd <- matrix(summary.chain$statistics,NCOL(x),4)[mask,2]
    x <- x[, mask, drop = FALSE]
    plot.priors <- FALSE
    if(!is.null(priors) && is.list(priors) &&  NCOL(priors[[1]])==NCOL(x) && NCOL(priors[[2]])==NCOL(x)) {
        priors[[1]] <- priors[[1]][, mask, drop = FALSE]
        priors[[2]] <- priors[[2]][, mask, drop = FALSE]
        plot.priors <- TRUE }
    
    data <- split(x, col(x), drop = FALSE)
    names(data) <- unlist(lapply(attributes(x)$dimnames[[2]][mask],
                                 function(z) chartr(":","_",z)))
    if(!is.na(sel <- match("(Intercept)",names(data))))
        names(data)[[sel]] <- "Intercept"
    xnam <- paste("~", names(data)[1])
    if(colx > 1) {
        for (j in 2:colx)
            xnam <- paste(xnam, names(data)[j], sep = "+") 
    }
    plotting <- as.formula(xnam)

    if (normal.density) {
        if (colx != length(arg.mean) | colx != length(arg.sd)) {
            warning("mean and sd of wrong dimension")
            print(hist)
            return }
        ynorm <- matrix()
        xnorm <- matrix(0,colx,2)
        for (j in 1:colx) {
            ynorm[j] <- dnorm(arg.mean[j], arg.mean[j], arg.sd[j])
            xnorm[j,] <- cbind(qnorm(0.00015,arg.mean[j],arg.sd[j]),
                               qnorm(0.00015,arg.mean[j],arg.sd[j],lower = FALSE))
        }
        ylim <- hist$y.limits
        xlim <- hist$x.limits
        for (j in 1:colx) {
            if (ylim[[j]][2] < ynorm[j])
                ylim[[j]][2] <- ynorm[j]
            if (xlim[[j]][1] > xnorm[j,1])
                xlim[[j]][1] = xnorm[j,1]
            if (xlim[[j]][2] < xnorm[j,2])
                xlim[[j]][2] = xnorm[j,2]
        }
        hist$x.limits <- xlim
        hist$y.limits <- ylim
      
        hist <- histogram( x = plotting, data = data, xlab = xlab, type = 
                          "density",
                          prepanel = function(x) list(xlim = range(x)),
                          panel = function(x, ...) {
                              panel.histogram(x, col = col,...)
                              panel.mathdensity(dmath = dnorm, col = density.col,
                                                lty = density.lty,
                                                args = list(mean = arg.mean[panel.number()],
                                                sd = arg.sd[panel.number()]))
                          },
                          layout = c(2,ceiling(colx/2) ), breaks = NULL,
                          nint = nint, scales = list(relation = "free",
                                       alternating = 3, col = "black"), ...)
        hist$x.limits <- xlim
        hist$y.limits <- ylim
    }
    else
    {
        if(kernel.est.density) {
            hist <- histogram( plotting, data = data,
                              xlab = xlab, type = "density",
                              panel = function(x, ...) {
                                  panel.histogram(x, ...)
                                  panel.lines(density(x)$x,density(x)$y,col = "red")
                              } )}
        else
            hist <- histogram( x = plotting, data = data, xlab = xlab, type = 
                              "density", 
                              prepanel = function(x) list(xlim = range(x)),
                              panel = function(x, ...) { panel.histogram(x, col = col,...) },
                              layout = c(2,ceiling(colx/2) ), breaks = NULL,
                              nint = nint, scales = list(relation = "free",
                                           alternating = 3, col = "black"), ...)

    }
    print(hist)
    if(plot.priors) {
        panel.locs <- trellis.currentLayout()
        i <- 1

        for (row in 1:nrow(panel.locs))
            for (column in 1:ncol(panel.locs))
                if (panel.locs[row, column] > 0)
                {
                    trellis.focus("panel", row = row, column = column,
                                  highlight = FALSE)
                    panel.lines(x = priors[[1]][,i], y = priors[[2]][,i], col = 'red')
                    trellis.unfocus()
                    i <- i + 1
                }

    }
}

tune.mcmc <- function( x, ... )
{
    attributes(x)$tune
}


