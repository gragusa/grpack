##################################################################
## utils.R /
## Author: Giuseppe Ragusa 
## Time-stamp: "2012-03-25 17:05:51 gragusa" 
##
## Description: Utils function for grpack
##################################################################
##' @export
rndCGM <- function( alpha = 0, beta = 1, sdalpha = 1, sdepsilon = 1,
                   nclusters = 4, nsubjects = 30, errortype = "homoskedastic",
                   errordist = 1, 
                   constantreg = FALSE, nins = 0, rho = 0.9, r2 = 0.1,
                   nuisance = 0, sdnuisance = 1)
{
    s  <- 1:nclusters
    S  <- nclusters
    Ns <- nsubjects
    if(errordist == 1)
        ag <- rep(rnorm(S, mean = 0, sd = sdalpha), each = Ns)
    else
        ag <- rep(sqrt(sdalpha)*(rchisq(S, df=1)-1)/sqrt(2), each = Ns)
    Xg <- rep(rnorm(S, mean = 0, sd = sdalpha), each = Ns)
    
    ## Generate Instruments
    
    if(!constantreg) {
        Xig <- rnorm(S*Ns, mean = 0, sd = sdalpha)
        scale <- 1
    }
    else {
        Xig <- rep(rnorm(S, mean = 0, sd = sdalpha), each = Ns)
        scale <- sqrt(2)
    }
    
    if (errortype == "homoskedastic") {
        eig <- rnorm(S*Ns,sd = sdepsilon)
        uig <- rnorm(S*Ns,sd = sdepsilon)
    } else {
        eig <- rnorm(S*Ns,sd = sqrt(9*(Xg+Xig)^2))
        uig <- rnorm(S*Ns,sd = sdepsilon)
    }
    
    if(nins > 0) {
        gamma = sqrt(r2/(nins*(1-r2)))
        errc <- cbind(eig,uig)
        varerrc <- cbind(c(1,rho), c(rho,1))
        nerr <- errc%*%chol(varerrc)
        eig <- nerr[,1]
        uig <- nerr[,2]
        if(!constantreg)
            Zig <- matrix(rnorm(nins*S*Ns, mean = 0, sd = sdalpha),S*Ns,nins)
        else
            Zig <- apply(as.matrix(1:ninst), 1,
                         function(x) rep(rnorm(S, mean = 0, sd = sdalpha), each = Ns))
        if(errordist == 1)
          nug <- rep(rnorm(S, mean = 0, sd = sdalpha), each = Ns)
        else
          nug <- rep(sqrt(sdalpha)*(rchisq(S, df=1)-1)/sqrt(2), each = Ns)    
        
        Xig <- Zig%*%rep(gamma,nins) + uig + nug
    }
    
    yis <- alpha+beta*((Xg+Xig)/scale)+eig+ag
    if(nins)
        data.frame(y = yis, X = (Xig+Xg)/scale, Zig = Zig, cluster = rep(1:S,each = Ns))
    else
        data.frame(y = yis, X = (Xig+Xg)/scale, cluster = rep(1:S,each = Ns))
}

rndCGM2 <- function( sdalpha = 1, k = 5, nclusters = 5, nsubjects = 50, 
                    gamma = c(1,.025,1,0) )
{
    k <- k-2
    s  <- 1:nclusters
    S  <- nclusters
    Ns <- nsubjects
    ag <- rep(rnorm(S, mean = 0, sd = sdalpha), each = Ns)    
    U <- matrix(rnorm(S*Ns*k), S*Ns, k)
    Xg <- rep(rnorm(S)*runif(S, min = 1, max =3)/2, each = Ns)
    ##     x1 <- rnorm(S*Ns, mean = 0.5, sd = 1.2)
    ##     x2 <- rnorm(S*Ns, mean = -0.5, sd = .7)
    ##     V <- .5*x1+.5*x2
    V <- rnormmix(S*Ns, lambda = c(.5,.5), mu = c(.5, -.5), sigma = c(1.2, 0.7))
    Z <- runif(S*Ns, min = 1, max = 3)
    X <- cbind(1,Xg,U*c(Z)/2, deparse.level = 0)
    Q <- (apply(X, 1, function(x) sum(x^2))-1-4.333)
    Y <- (gamma[1]*Q+gamma[2]*Q*c(V)*.25+gamma[3]*V+gamma[4]*ag)/((4+.2)*gamma[1]+(1-gamma[1])*sqrt(2))
    data.frame(Y = Y, X = X[,-1], cluster = rep(1:S,each = Ns))
}


trim <- function(x, q = 0.01) {
  qq <- quantile(x, p=c(q, 1-q))
  x[x>qq[1] & x<qq[2]]
}

capply <- function(str, ff) 
  sapply(lapply(strsplit(str, NULL), ff), paste, collapse="") 

cap <- function(char) {
  ## change lower letters to upper, others leave unchanged
  if (any(ind <- letters==char)) LETTERS[ind]    else char 
}

capitalize <- function(str) { # vector of words
   ff <- function(x) paste(lapply(unlist(strsplit(x, NULL)),cap),collapse="")
   capply(str,ff) 
}

lower <- function(char) {
    ## change upper letters to lower, others leave unchanged
    if (any(ind <- LETTERS==char)) letters[ind]    else char 
}

lowerize <- function(str) {
    ff <- function(x) paste(lapply(unlist(strsplit(x, NULL)),lower),collapse="")
    capply(str,ff) 
}

"CapLeading" <- function(str) {
    ff <- function(x) {r <- x; r[1]<-cap(x[1]); r}
    capply(str,ff) 
}

model.cluster <- function (x) 
  x$"(cluster)"
