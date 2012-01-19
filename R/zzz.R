##################################################################
## zzz.R /
## Author: Giuseppe Ragusa 
## Time-stamp: "2011-11-02 19:01:41 gragusa" 
##
## Description:
##################################################################

.First.lib <-
    function(lib, pkg) {
        library.dynam("grpack", pkg, lib)
        require(stats, quietly = TRUE)    
    }
