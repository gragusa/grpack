##################################################################
## hacks.R /
## Author: Giuseppe Ragusa 
## Time-stamp: "2010-03-09 16:02:50 gragusa hacks.R" 
##
## Description: Collect useful hacks
##################################################################


Title <- function(...)
{
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(mfrow = c(1, 1), xpd = NA)
    plot.window(c(0,1),c(0,1))
    title(...)
}

