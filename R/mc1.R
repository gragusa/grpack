##################################################################
## mc1.R /
## Author: Giuseppe Ragusa 
## Time-stamp: "2008-11-07 16:30:07 gragusa" 
##
## Description: 
##################################################################

mc1 <- function(simul = 100, nclus = 10, true = 1, n0 = 1, bootsim = 399, geweights = 'radamacher',
                wbweights = 'radamacher', errortype = 'homoskedastic', errordist = 1)
{
    out <- matrix(0, simul, 5)
    for(j in 1:simul)
    {
        dat <- rndCGM(beta = true, ncluster = nclus, errortype = errortype, errordist = errordist)
        full <- reg(y~X, data = dat, cluster = cluster)
        ge <- geboot.reg(full, sim = bootsim, wbweights = 'exp')
        wb <- wildboot.reg(full, sim = bootsim, null = list(X=n0), wbweight = 'radamacher')
        sge <- summary(ge, null = list(X=n0))
        swb <- summary(wb, null = list(X=n0))
        out[j, ] <- c(sge$pvalue, sge$pvalue.1, swb$pvalue, swb$pvalue.1, swb$normal.pvalue)
    }
    out
}
