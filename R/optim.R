##################################################################
## optim.R /
## Author: Giuseppe Ragusa 
## Time-stamp: "2012-11-03 19:08:12 gragusa" 
##
## Description: optimization routines:
## bfgs -- smooth function minimization with the BFGS algorithm
## based on Andrew Clausen <clausen@econ.upenn.edu> code
## coslve -- equation solver based on Chris Sims
##################################################################
interpolate.min <- function(x0, x1, f0, f1, g0)
{
	dx <- x0 - x1
	dx2 <- x0^2 - x1^2

	a <- (f0 - f1 - g0*dx) / (dx2 - 2*x0*dx)
	b <- g0 - 2*a*x0

	x <- -2*a / b
	if (is.na(x) || (x < min(x0, x1)) || (x > max(x0, x1)))
		return ((x0 + x1) / 2)
	x
}

linesearch <- function(phi, phi_, alpha1, alpha.max,
		       ftol=0.0001, gtol=0.9, stepsize=3)
{
    ## Implements the More-Thuente linesearch algorithm.
    ## The notation is taken from
    ## Nocedal and Wright (2006).
    ## This function returns an approximate solution to
    ##
    ##	argmin_{alpha \in [0, alpha.max]} phi(alpha).
    
    ## the Wolfe conditions are: (we want both to be TRUE)
    armijo <- function(alpha)
        phi(alpha) < phi(0) + ftol * alpha * phi_(0)
    curvature <- function(alpha)
        abs(phi_(alpha)) <= gtol * abs(phi_(0))
    
    zoom <- function(alpha.lo, alpha.hi) {
        for (i in 1:30) {
            alpha <- interpolate.min(
                                     alpha.lo, alpha.hi, phi(alpha.lo),
                                     phi(alpha.hi), phi_(alpha.lo))
            
            if (!armijo(alpha) || (phi(alpha) >= phi(alpha.lo))) {
                alpha.hi <- alpha
            } else {
                if (curvature(alpha))
                    return(alpha)
                                        # started going uphill?
                if (phi_(alpha) * (alpha.hi - alpha.lo) >= 0)
                    alpha.hi <- alpha.lo
                alpha.lo <- alpha
            }
        }
        0	# not enough progress; give up.
    }
    
    stopifnot(phi_(0) < 0)
    
    alpha_ <- 0
    alpha <- ifelse(alpha1 >= alpha.max, alpha.max/2, alpha1)
    for (i in 1:100) {
        if (i > 1 && (phi(alpha) >= phi(alpha_)))
            return(zoom(alpha_, alpha))
        if (!armijo(alpha))
            return(zoom(alpha_, alpha))
        if (curvature(alpha))
            return(alpha)
        if (phi_(alpha) >= 0)
            return(zoom(alpha, alpha_))
        alpha_ <- alpha
        alpha <- min((alpha + alpha.max) / 2, alpha * stepsize)
    }
}

norma <- function(x) max(abs(x))

# collects statistics on how many times f is called.
call.counter <- function(f, name, environment)
{
	assign(name, 0, envir=environment)
	function(x)
	{
		assign(name, 1 + get(name, envir=environment),
		       envir=environment)
		f(x)
	}
}

## Returns a function wrapper of f that caches old values.
## (i.e. memorization)
function.cache <- function(f)
{
    cache.env <- new.env()
    cache <- list()
    cache$n <- 0
    cache$params <- list()
    cache$vals <- list()
    assign("cache", cache, envir=cache.env)
    
    function(x)
    {
        cache <- get("cache", envir=cache.env)
        if (cache$n > 0) {
            for (i in cache$n:max(1, (cache$n-100))) {
                if (all(x == cache$params[[i]]))
                    return(cache$vals[[i]])
            }
        }
        cache$n <- cache$n + 1
        cache$params[[cache$n]] <- x
        cache$vals[[cache$n]] <- f(x)
        assign("cache", cache, envir=cache.env)
        cache$vals[[cache$n]]
    }
}

# Wrapper for numericDeriv that returns the derivative function of g.
# (numericDeriv only evalutes the derivative.)
numericDeriv_ <- function(f, grad=FALSE)
{
    rho <- new.env()
    assign("f", f, envir=rho)
    function(x)
    {
        assign("x", x, envir=rho)
        result <- attr(numericDeriv(quote(f(x)), "x", rho), "gradient")
        if (grad)
            result <- as.vector(result)
        result
    }
}

## We require each coordinate of x satisfy
#
#	x in [min.x, max.x].
#
## This function returns max.lambda such that
## for all lambda in [0, max.lambda],
##
##	(x + s * max.lambda) in [min.x, max.x].
calc.constraint <- function(x, s, min.x, max.x)
{
    stopifnot((x >= min.x) && (x <= max.x))
    
    min.constraints <- (x - min.x) / (-s)
    min.constraints <- subset(min.constraints, min.constraints > 0)
    
    max.constraints <- (max.x - x) / s
    max.constraints <- subset(max.constraints, max.constraints > 0)
    
    constraints <- c(min.constraints, max.constraints)
    if (all(is.infinite(constraints)))
        return(Inf)
    min(constraints)
}

## Implements the Broyden-Fletcher-Goldfarb-Shanno algorithm
## for function minimization.  That is, it attempts to find
##
##	argmin_{x s.t. min.x < x < max.x coord-wise} f(x)
##
## Apart from returning the maximizer x, it also returns the
## function value f(x)
## the gradient f'(x), and the inverse hessian f''(x)^{-1}.

##' Implements the Broyden-Fletcher-Goldfarb-Shanno algorithm for function
##' minimization.
##' 
##' The function attempts to find
##'
##' argmin_{x s.t. min.x < x < max.x coord-wise} f(x)
##'
##' 
##' @title bfgs
##' @param x0 Tthe initial solution guess
##' @param f_ The function to be minimized
##' @param g_ The gradient of f_
##' @param min.x lower bounds
##' @param max.x upper bounds
##' @param prec 
##' @param verbose if TRUE diplay iteration nformation
##' @return A list
##'
##' Apart from returning the maximizer x, it also returns the
##' function value f(x) the gradient f'(x), and the inverse hessian
##' f''(x)^{-1}.

##' @author Giuseppe Ragusa
bfgs <- function(x0, f_, g_,
		 min.x=rep(-Inf, length(x0)),
		 max.x=rep(Inf, length(x0)),
		 prec=0.00001, verbose=FALSE)
{
    count.env <- new.env()
    f <- function.cache(call.counter(function(x) f_(x), "f", count.env))
    
    if (missing(g_))
        g_ <- numericDeriv_(f, grad=TRUE)

    g <- function.cache(call.counter(function(x) g_(x), "g", count.env))
    
    
    x <- x0
    I <- diag(rep(1, length(x)))
    H <- I
    convergence <- 0
    iter <- iterf <- 0
    if(class(try(g(x), silent = TRUE))=='try-error')
        while (class(try(g(x), silent = TRUE))=='try-error' & iterf<200) {
            iterf <- iterf+1
            x <- x+runif(length(x))
        }
    while (norma(g(x)) > prec) {
        iter <- iter + 1
        if (verbose)
            cat(c("\niter", iter, f(x), x, "\n"))
        
        ## minimize in the direction of p
        p <- as.vector(- H %*% g(x))
        phi <- function(alpha) f(x + alpha * p)
        phi_ <- function(alpha) as.numeric(g(x + alpha * p) %*% p)
        max.alpha <- calc.constraint(x, p, min.x, max.x)
        alpha <- linesearch(phi, phi_, 1, max.alpha)
        if (alpha == 0) {
            ##warning("Lost precision")
            convergence <- 1
            break
        }
        
        x_ <- x + alpha * p
        s <- x_ - x
        y <- g(x_) - g(x)
        
        rho <- 1 / as.numeric(t(y) %*% s)
        J <- I - s %*% t(y) * rho
        K <- I - y %*% t(s) * rho
        H <- J %*% H %*% K + s %*% t(s) * rho
        
        x <- x_
    }
    
## 	list(par=x, value=f(x), grad=g(x), inv.hessian=H,
## 	     counts=c(`function`=get("f", envir=count.env),
## 		      gradient=get("g", envir=count.env)),
##              convergence = convergence)
    list(par=x, value=f(x), grad=g(x), inv.hessian=H,
         convergence = convergence, iter = iter)
}

bfgs.2 <- function(x0, f_, g_, mask,
                   min.x=rep(-Inf, length(x0)),
                   max.x=rep(Inf, length(x0)),
                   prec=0.00001, verbose=FALSE)
{
    xf <- numeric()
    x <- x0[!mask$pres]
    f <- function(x) {
      xf <- numeric(length=length(x0))
      xf[mask$pres] <- mask$pnull
      xf[!mask$pres] <- x
      f_(xf)
    }
    g <- function(x) {
      xf <- numeric(length=length(x0))
      xf[mask$pres] <- mask$pnull
      xf[!mask$pres] <- x
      g_(xf)[,!mask$pres]
    }
    
    I <- diag(rep(1, length(x)))
    H <- I
    iter <- 0
    
    while (norma(g(x)) > prec) {
      iter <- iter + 1
      if (verbose)
        cat(c("\niter", iter, f(x), x, "\n"))
      
      ## minimize in the direction of p
      p <- as.vector(- H %*% gr)
      phi <- function(alpha) f(x + alpha * p)
      phi_ <- function(alpha) as.numeric(g(x + alpha * p) %*% p)
      max.alpha <- calc.constraint(x, p, min.x, max.x)
      alpha <- linesearch(phi, phi_, 1, max.alpha)
      if (alpha == 0) {
        warning("Lost precision")
        break
      }
      
      x_ <- x + alpha * p
      s <- x_ - x
      y <- g(x_) - gr
      
      rho <- 1 / as.numeric(t(y) %*% s)
      J <- I - s %*% t(y) * rho
      K <- I - y %*% t(s) * rho
      H <- J %*% H %*% K + s %*% t(s) * rho
      
      x <- x_
      }
    
    list(par=x, value=f(x), grad=g(x), inv.hessian=H)
  }

csolve <- function(FUN,FUN1,x,...,gradfun = NULL,
                   crit = 1e-6,
                   itmax = 10,
                   verbose = TRUE,
                   alpha = 1e-3,
                   delta = 1e-6,
                   long = FALSE)
{

  EPS <- .Machine$double.eps
  if (is.null(dim(x))) x <- matrix(x,length(x),1)
  nv <- dim(x)[1]
  vdn <- dimnames(x)
  tvec <- delta*diag(nv)
  done <- FALSE
  f0 <- FUN(x,...)
  grad <- f0[[2]]
  f0 <- f0[[1]]
  af0 <- sum(abs(f0))
  af00 <- af0
  itct <- 0
  while(!done) {
    if((itct%%2) == 1  && af00-af0 < crit*max(1,af0) && itct > 3) {
      randomize <- TRUE
    } else {
      if(is.finite(as.matrix(grad)) && sum(abs(grad)) > 4*nv^2*EPS) {
        svdg <- svd(grad)
        svdd <- svdg$d
        if(!(min(svdd) > 0) || max(svdd)/min(svdd) > 100*EPS){
          svdd <- pmax(svdd,max(svdd)*1e-13)
          grad <- svdg$u%*% diag(svdd,ncol = length(svdd)) %*% t(svdg$v)
        }
        dx0 <- -solve(grad,f0)
        randomize <- FALSE
      } else {
        ##      if(verbose){cat("gradient imaginary or infinite or zero\n")}
        randomize <- TRUE
      }
    }
    if(randomize) {
      ##      if(verbose){cat("Random Search\n")}
      dx0 <- sqrt(sum(x*x))/matrix(rnorm(length(x)),length(x),1)
    }
    lambda <- 1
    lambdamin <- 1
    fmin <- f0
    xmin <- x
    afmin <- af0
    dxSize <- sqrt(sum(abs(dx0)^2))
    factor <- .6
    shrink <- TRUE
    subDone <- FALSE
    while(!subDone) {
      dx <- lambda*dx0
      f <- FUN1(x+dx)
      af <- sum(abs(f))
      if(af < afmin) {
        afmin <- af
        fmin <- f
        lambdamin <- lambda
        xmin <- x+dx
      }
      if( ((lambda > 0) && (af0-af < alpha*lambda*af0)) || ((lambda < 0) && (af0-af < 0) )) {
        if(!shrink) {
          factor <- factor^.6
          shrink <- TRUE
        }
        if(abs(lambda*(1-factor))*dxSize > .1*delta) {
          lambda <- factor*lambda
        } else {
          if( (lambda > 0) && (factor == .6) ) { #i.e., we've only been shrinking
            lambda <- -.3
          } else {
            subDone <- TRUE
            if(lambda > 0) {
              if(factor == .6) {
                rc <- 2
              } else {
                rc <- 1
              }
            } else {
              rc <- 3
            }
          }
        }
      } else {
        if( (lambda > 0) && (af-af0 > (1-alpha)*lambda*af0) ) {
          if(shrink) {
            factor <- factor^.6
            shrink <- FALSE
          }
          lambda <- lambda/factor
        } else {                        # good value found
          subDone <- TRUE
          rc <- 0
        }
      }
    }                                   # while ~subDone
    itct <- itct+1
    x <- xmin
    f0 <- fmin
    af00 <- af0
    af0 <- afmin
    if(itct >= itmax) {
      done <- TRUE
      rc <- 4
    } else {
      if(af0 < crit) {
        done <- TRUE
        rc <- 0
      }
    }
  }
  
  return(list(par = xmin,value = 0, gradient = f0, hessian = grad, convergence = rc))
}

csolve <- function(FUN,FUN1,x,...,gradfun = NULL,crit = 1e-6,itmax = 10,
                   verbose = TRUE,alpha = 1e-3,delta = 1e-6,long = FALSE)
{

  EPS <- .Machine$double.eps
  if (is.null(dim(x))) x <- matrix(x,length(x),1)
  nv <- dim(x)[1]
  vdn <- dimnames(x)
  tvec <- delta*diag(nv)
  done <- FALSE
  f0 <- FUN(x,...)
  grad <- f0[[2]]
  f0 <- f0[[1]]
  af0 <- sum(abs(f0))
  af00 <- af0
  itct <- 0
  while(!done) {
    if((itct%%2) == 1  && af00-af0 < crit*max(1,af0) && itct > 3) {
      randomize <- TRUE
    } else {
      if(is.finite(as.matrix(grad)) && sum(abs(grad)) > 4*nv^2*EPS) {
        svdg <- svd(grad)
        svdd <- svdg$d
        if(!(min(svdd) > 0) || max(svdd)/min(svdd) > 100*EPS){
          svdd <- pmax(svdd,max(svdd)*1e-13)
          grad <- svdg$u%*% diag(svdd,ncol = length(svdd)) %*% t(svdg$v)
        }
        dx0 <- -solve(grad,f0)
        randomize <- FALSE
      } else {
        ##      if(verbose){cat("gradient imaginary or infinite or zero\n")}
        randomize <- TRUE
      }
    }
    if(randomize) {
      ##      if(verbose){cat("Random Search\n")}
      dx0 <- sqrt(sum(x*x))/matrix(rnorm(length(x)),length(x),1)
    }
    lambda <- 1
    lambdamin <- 1
    fmin <- f0
    xmin <- x
    afmin <- af0
    dxSize <- sqrt(sum(abs(dx0)^2))
    factor <- .6
    shrink <- TRUE
    subDone <- FALSE
    while(!subDone) {
      dx <- lambda*dx0
      f <- FUN1(x+dx)
      af <- sum(abs(f))
      if(af < afmin) {
        afmin <- af
        fmin <- f
        lambdamin <- lambda
        xmin <- x+dx
      }
      if( ((lambda > 0) && (af0-af < alpha*lambda*af0)) || ((lambda < 0) && (af0-af < 0) )) {
        if(!shrink) {
          factor <- factor^.6
          shrink <- TRUE
        }
        if(abs(lambda*(1-factor))*dxSize > .1*delta) {
          lambda <- factor*lambda
        } else {
          if( (lambda > 0) && (factor == .6) ) { #i.e., we've only been shrinking
            lambda <- -.3
          } else {
            subDone <- TRUE
            if(lambda > 0) {
              if(factor == .6) {
                rc <- 2
              } else {
                rc <- 1
              }
            } else {
              rc <- 3
            }
          }
        }
      } else {
        if( (lambda > 0) && (af-af0 > (1-alpha)*lambda*af0) ) {
          if(shrink) {
            factor <- factor^.6
            shrink <- FALSE
          }
          lambda <- lambda/factor
        } else {                        # good value found
          subDone <- TRUE
          rc <- 0
        }
      }
    }                                   # while ~subDone
    itct <- itct+1
    x <- xmin
    f0 <- fmin
    af00 <- af0
    af0 <- afmin
    if(itct >= itmax) {
      done <- TRUE
      rc <- 4
    } else {
      if(af0 < crit) {
        done <- TRUE
        rc <- 0
      }
    }
  }
  return(list(par = xmin,value = 0, gradient = f0, hessian = grad, convergence = rc))
}
