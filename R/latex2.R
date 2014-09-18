##' Exactly like latex of \code{Hmisc}, but fix a bug when using
##' \code{ctable=TRUE} and add three arguments: \code{do.begin},
##' \code{do.end}, \code{tight}. 
##'
##' Details here
##' 
##' @title Latex tables.
##' @param do.begin Should the latex code at the beginning of table be printed?
##' @param do.end Should the latex code at the end of table be printed?
##' @param tight when using dcolumn tighten columns width.
##' @param rowcolors This command is useful for inserting color command. For
##' instance, if \code{rowcolors=rowcolors{1}{green}{pink}}, then the rows are 
##' subsequently coloured in green and in pink. 
##' @return See \code{\link{latex}} of package \code{Hmisc}.
##' @author Giuseppe Ragusa
##' 
##' @export
latex2 <-
  function(object,
           title=first.word(deparse(substitute(object))),
           file=paste(title, ".tex", sep=""),
           append=FALSE, do.begin = TRUE,
           do.end=TRUE, label=title,
           rowlabel=title, rowlabel.just="l",
           cgroup=NULL, n.cgroup=NULL,
           rgroup=NULL, n.rgroup=NULL,
           cgroupTexCmd="bfseries",
           rgroupTexCmd="bfseries",
           rownamesTexCmd=NULL, 
           colnamesTexCmd=NULL,
           cellTexCmds=NULL,
           rowname, cgroup.just=rep("c",sum(n.cgroup)),
           colheads=NULL,
           extracolheads=NULL, extracolsize='scriptsize',
           dcolumn=FALSE, tight=FALSE, numeric.dollar=!dcolumn, cdot=FALSE,
           longtable=FALSE, draft.longtable=TRUE, ctable=FALSE, booktabs=FALSE,
           table.env=TRUE, here=FALSE, lines.page=40,
           caption=NULL, caption.lot=NULL, caption.loc=c('top','bottom'),
           double.slash=FALSE,
           vbar=FALSE, collabel.just=rep("c",nc), na.blank=TRUE,
           insert.bottom=NULL, first.hline.double=!(booktabs | ctable),
           where='!tbp', size=NULL,
           center=c('center','centering','none'),
           landscape=FALSE,
           multicol=TRUE, ## to remove multicolumn if no need  SSJ 17nov03
           math.row.names=FALSE, math.col.names=FALSE, rowcolors=NULL,
           ...)      ## center MJ 08sep03
{
  center <- match.arg(center)
  caption.loc <- match.arg(caption.loc)
  cx <- format.df(object, dcolumn=dcolumn, na.blank=na.blank,
                  numeric.dollar=numeric.dollar, cdot=cdot,
                  math.row.names=math.row.names, math.col.names=math.col.names,
                  double.slash=double.slash, ...)
  ## removed check.names=FALSE from above 23jan03
  if (missing(rowname))
    rowname <- dimnames(cx)[[1]]
  
  if (is.null(colheads))
    colheads <- dimnames(cx)[[2]]

  col.just <- attr(cx,"col.just")
  nc <- ncol(cx)
  nr <- nrow(cx)

  if (length(cgroup)) {
    k <- length(cgroup)
    if(!length(n.cgroup))
      n.cgroup <- rep(nc/k, k)
    
    if(sum(n.cgroup)!=nc)
      stop("sum of n.cgroup must equal number of columns")
    
    if(length(n.cgroup)!=length(cgroup))
      stop("cgroup and n.cgroup must have same lengths")
  }

  if(!length(rowname))
    rgroup <- NULL
  
  if(!length(n.rgroup) && length(rgroup))
    n.rgroup <- rep(nr/length(rgroup), length(rgroup))
  
  if(length(n.rgroup) && sum(n.rgroup)!=nr)
    stop("sum of n.rgroup must equal number of rows in object")
  
  if(length(rgroup) && length(n.rgroup) && (length(rgroup)!=length(n.rgroup)))
    stop("lengths of rgroup and n.rgroup must match")
  
  if (length(rgroup) && rowlabel.just=="l")
    rowname <- paste("~~",rowname,sep="")

  sl <- ifelse(double.slash, "\\\\", "\\")
  if(ctable) {
    eol <- paste(sl, 'NN\n', sep='')
    eog <- ""
  } else if(longtable && length(n.rgroup)) {
    eol <- paste(sl,"tabularnewline*\n",sep='')
    eog <- paste(sl, "tabularnewline\n", sep='')
  } else {
    eol <- paste(sl,"tabularnewline\n",sep='')
    eog <- paste(sl, "tabularnewline\n", sep='')      
  }
  
  if(booktabs) {  # 27may02
    toprule    <- paste(sl,"toprule\n",sep="")
    midrule    <- paste(sl,"midrule\n",sep="")
    bottomrule <- paste(sl,"bottomrule\n",sep="")
  } else if(ctable) {   ## 13dec02
    toprule    <- paste(sl, 'FL\n', sep='')
    midrule    <- paste(sl, 'ML\n', sep='')
    bottomrule <- paste(sl, 'LL\n', sep='')
  } else {
    toprule <-
      if(first.hline.double)
        paste(sl,"hline",sl,"hline\n",sep="")
      else
        paste(sl,"hline\n",sep="")
    
    midrule <- bottomrule <- paste(sl,"hline\n",sep="")
  }


  ## ################ CELL AND ROWNAMES FORMATS ###################
  ## If no formats are specified for the rownames and cells there is
  ## nothing to do. If only one is specified then the other must
  ## faked. But rownamesTexCmd should only be faked if rownames is
  ## not NULL.

  ## Check to make sure the dimensions of the cell formats
  ## match the dimensions of the object to be formatted.
  if (!is.null(cellTexCmds) &
      !(all(dim(cx) == dim(cellTexCmds)) &
        length(dim(cx)) == length(dim(cellTexCmds)))) {
    msg <- "The dimensions of cellTexCmds must be:"
    msg1 <- paste(dim(cx), collapse=" x ")
    msg <- paste(msg, msg1)
    msg <- paste(msg, ", but you gave me: ")
    msg1 <- paste(dim(cellTexCmds), collapse=" x ")
    msg <- paste(msg, msg1, sep="")
    stop(msg)
  }
  
  ## If there are column groups, add a blank column
  ## of formats between the groups.
  if (length(cgroup) & !is.null(cellTexCmds)) {
    my.index <- split(1:NCOL(cellTexCmds), rep(cumsum(n.cgroup), times=n.cgroup))
    new.index <- NULL
    new.col <- dim(cx)[2] + 1
    for (i in my.index)
      new.index <- c(new.index, i, new.col)
    
    new.index <- new.index[-length(new.index)]
    cellTexCmds <- cbind(cellTexCmds, "")[, new.index]
  }

  if (!is.null(cellTexCmds) | !is.null(rownamesTexCmd)) {
    ## LaTeX commands have been specified for either the rownames or
    ## the cells.
    ## Fake rownamesTexCmd if it is NULL and if rowname exists.
    if (is.null(rownamesTexCmd) & !is.null(rowname))
      rownamesTexCmd <- rep("", nr)
    
    ## Fake cellTexCmds if it is NULL.
    if (is.null(cellTexCmds)) {
      cellTexCmds <- rep("", dim(cx)[1] * dim(cx)[2])
      dim(cellTexCmds) <- dim(cx)
    }
    
    ## Create a combined rowname and cell format object.
    rcellTexCmds <- cbind(rownamesTexCmd, cellTexCmds)
    thisDim <- dim(rcellTexCmds)
    ## Prefix the latex commands with slashes.
    rcellTexCmds <- paste(sl, rcellTexCmds, sep="")
    ## Remove slashes from elements where no format was specified.
    rcellTexCmds[rcellTexCmds == sl] <- ""
    ## Restore the dimensions of the matrix (paste loses them).
    dim(rcellTexCmds) <- thisDim
  } else {
    rcellTexCmds <- NULL
  }

  ## ############## END OF CELL AND ROWNAMES FORMATS ###############
  
  
  ##if (!vbar && length(cgroup)) {
  if (length(cgroup)) {
    last.col <- cumsum(n.cgroup)
    first.col <- c(1, 1+last.col[-length(last.col)])
    cgroup.cols <- cbind(first.col,last.col)
    col.subs <- split(seq(length.out=nc), rep.int(seq_along(n.cgroup), times=n.cgroup))
    
    cxi <- list()
    for (i in seq(along=col.subs))
      cxi[[i]] <- cx[,col.subs[[i]],drop=FALSE]
    
    cxx <- cxi[[1]]
    col.justxx <- col.just[col.subs[[1]]]
    collabel.justxx <- collabel.just[col.subs[[1]]]
    colheadsxx <-  colheads[col.subs[[1]]]
    extracolheadsxx <- extracolheads[col.subs[[1]]]

    cgroupxx <- cgroup[1]
    n.cgroupxx <- n.cgroup[1]
    for (i in seq(along=col.subs)[-1]) {
      cxx <- cbind(cxx, "", cxi[[i]])  # was ""="" 23Feb01 "=" 2Apr02
      col.justxx <- c(col.justxx, "c", col.just[col.subs[[i]]])
      collabel.justxx <- c(collabel.justxx, "c",
                           collabel.just[col.subs[[i]]])
      cgroupxx <- c(cgroupxx, "", cgroup[i])
      n.cgroupxx <- c(n.cgroupxx, 1, n.cgroup[i])
      colheadsxx <- c(colheadsxx, "", colheads[col.subs[[i]]])
      if(!is.null(extracolheads)) {
        extracolheadsxx <- c(extracolheadsxx, "",
                             extracolheads[col.subs[[i]]])
      }
    }
    
    cgroup.colsxx <- cgroup.cols + 0:(nrow(cgroup.cols)-1)
    
    cx <- cxx
    col.just <- col.justxx
    collabel.just <- collabel.justxx
    n.cgroup <- n.cgroupxx
    cgroup.cols <- cgroup.colsxx[cgroup!="",,drop=FALSE]
    cgroup <- cgroupxx
    colheads <- colheadsxx
    extracolheads <- extracolheadsxx
    nc <- ncol(cx)
  }

  cline <- NULL
  if (length(rowname)) {
    cx <- cbind(rowname, cx)
    col.just <- c(rowlabel.just, col.just)

    if(length(extracolheads))
      extracolheads <- c('', extracolheads)  ## 16jun03
    
    collabel.just <- c(rowlabel.just, collabel.just)
    if (length(cgroup) == 0L)
      colheads <- c(rowlabel, colheads)
    else {
      colheads <- c('', colheads)
      cgroup <- c(rowlabel, cgroup)

      rlj <- ifelse(rowlabel.just=="l", "l", "c")
      cgroup.just <- c(rlj, cgroup.just)
      n.cgroup <- c(1, n.cgroup)
      cgroup.cols <- 1+cgroup.cols
      cline <- paste(sl, "cline{", cgroup.cols[,1],"-", cgroup.cols[,2], "}",
                     sep="", collapse=" ")
    }
    
    nc <- 1 + nc
  }

  vbar <- ifelse(vbar, "|", "")

  if(!append)
    cat("", file=file)	#start new file
  
  cat("%",deparse(sys.call()), "\n%\n", file=file, append=file!='')
  ## append= 19apr03 and other places
  ## Was as.character(as.name(match.call()))  15Sep00

  if(dcolumn) {
    decimal.point <- ifelse(cdot, paste(sl,"cdot",sep=""), ".")
    if(do.begin) {
      if(!tight) {
        cat(sl,"newcolumntype{.}{D{.}{",decimal.point,"}{-1}}\n",
            sep="", file=file, append=file!='')  # was newcolumn 26Feb02
      } else {
        cat(sl,"newcolumntype{.}{D{.}{",decimal.point,"}{1}}\n",
            sep="", file=file, append=file!='')  # was newcolumn 26Feb02
      }
    }
  }

  { # tabular.cols
    tabular.cols <- paste(vbar, col.just, sep="")
    if (!length(n.cgroup))
      tabular.cols <- c(tabular.cols, vbar)
    else {
      vv2 <- cumsum(n.cgroup)
      tabular.cols[vv2] <- paste(tabular.cols[vv2],vbar,sep="")
    }
    
    tabular.cols <- paste(tabular.cols, collapse="")
  }

  if(length(caption) && !ctable) {
    caption <- paste(sl,"caption",
                     if(length(caption.lot))
                       paste("[",caption.lot,"]",sep=""),
                     "{", caption,
                     if(!longtable)
                       paste(sl,"label{", label, "}",sep=""),
                     "}", sep="")
    
    table.env <- TRUE
  }

  if(ctable) {  ## 13dec02
    latex.begin <- paste(if(length(size))
                       paste('{',sl,size,sep=''),
                     paste(sl, "ctable[", sep=''),
                     if(length(caption) && caption.loc=='bottom')
                       'botcap,',
                     if(length(caption))
                       paste('caption={',caption,'},',sep=''),
                     if(length(caption.lot))
                       paste('cap={',caption.lot,'},',sep=''),
                     paste('label=',label,',',sep=''),
                     if(!landscape)
                       paste('pos=',where,',',sep=''),
                     if(landscape)
                       'sideways',
                     paste(']{',tabular.cols, '}',sep=''),
                     if(length(insert.bottom))
                       paste('{',
                             paste(sl,'tnote[]{',sedit(insert.bottom,'\\\\',' '),'}',
                                   sep='', collapse=''),
                             '}',
                             sep='')
                     else '{}',
                     ## tnote does not allow \\ in its argument
                     paste('{', toprule, sep=''),
                         sep='')
    
    latex.end <- paste('}',
                   if(length(size))
                     '}', sep='')
    
  } else if(!longtable) {
    latex.begin <- paste(if(table.env)
                         paste(sl, if (landscape) "begin{sidewaystable}" else "begin{table}",
                             if(here)
                               "[H]"
                             else
                               paste('[',where,']',sep=''),
                             "\n", sep=""),
                     if(length(size))
                       paste(sl,size,'\n',sep=''),
                     if(caption.loc=='top' && !missing(caption))
                       paste(caption, "\n"),              ## 3oct03
                     if(center == 'center')             ## MJ: 08sep03
                       paste(sl,"begin{center}\n", sep="")## MJ: 08sep03
                     else {
                       if (center == 'centering')  ## MJ: 08sep03
                         paste(sl,"centering\n", sep="")
                     }, ## MJ: 08sep03
                     paste(sl,"begin{tabular}{", tabular.cols, "}\n",
                           toprule, sep=""),
                     ## 11Jun95   12jan03 "}" was "}{" WHY!
                         sep='')
    
    latex.end <- paste(paste(sl,"end{tabular}\n", sep = ""),
                   if(center == 'center')  ## MJ: 08sep03
                     paste(sl,"end{center}\n", sep=""), ## MJ: 08sep03
                   if(caption.loc=='bottom' && !missing(caption))
                     paste(caption,'\n'),   # 3oct03
                   if(length(insert.bottom))
                     paste(insert.bottom, collapse='\\\\'),
                   if(table.env)
                     paste(sl, if(landscape) "end{sidewaystable)" else "end{table}\n", sep=""),
                   sep='')
  } else {
    latex.begin <- paste(paste(if (!draft.longtable)
                                 paste(sl,"let",sl,"LTmulticolumn=",sl,"multicolumn", sep=""),
                               paste(sl,"setlongtables",sep=""),
                                     if(landscape)
                                       paste(sl, "begin{landscape}",sep=""),
                               if(length(size))
                                 paste('{',sl,size,'\n',sep=''),
                               paste(sl,"begin{longtable}{", tabular.cols, "}",
                                     sep=""),
                               sep="\n"),
                         if(caption.loc=='top' && !missing(caption))
                           paste(caption, eog)
                         else 
                           '\n',
                         toprule, sep="")    #11Jun95
    
    latex.end <- paste(if(caption.loc=='bottom' && !missing(caption))
                         paste(caption, eog),  ## 3oct03
                       paste(sl,"end{longtable}\n", sep=""),
                       if(length(size))
                         '}\n',
                       if(landscape)
                         paste(sl,"end{landscape}\n",sep=""),
                       sep='')
  }
  if(do.begin)
      cat(latex.begin, file=file, append=file!='')

  if(length(cgroup)) {  # was !missing 5Oct00
    cvbar <- paste(cgroup.just, vbar, sep="")
    cvbar[1] <- paste(vbar, cvbar[1], sep="")
    cvbar[-length(cvbar)] <- paste(cvbar[-length(cvbar)], vbar, sep="")
    slmc <- paste(sl,"multicolumn{",sep="")
    ##labs <- paste(sl, "bf ", cgroup, sep="") 
    if (!is.null(cgroupTexCmd))
      labs <- paste(sl, cgroupTexCmd, " ", cgroup, sep="")
                                        # DRW 12apr05.
    
    if(multicol) ## SSJ 17nov03
      labs <- paste(slmc, n.cgroup, "}{", cvbar, "}{", labs, "}", sep="")

    cat(labs, file=file, sep="&\n", append=file!='')
    
    if (!length(cline)) {   # was is.length 2Apr02
      inr <- as.numeric(length(rowname))
      cline <- paste(sl,"cline{",1+inr,"-",nc,"}",sep="")
    }
    
    cat(eol, cline,"\n", sep="",file=file, append=file!='')
    ## eol was sl, sl  13dec02
  }


  { # column labels
    cvbar <- paste(collabel.just, vbar, sep="")
    cvbar[1] <- paste(vbar, cvbar[1], sep="")
    if (length(n.cgroup)) {
      vv2 <- cumsum(n.cgroup[-length(n.cgroup)])
      cvbar[vv2] <- paste(cvbar[vv2],vbar,sep="")
    }
    slmc1 <- paste(sl, "multicolumn{1}{", sep="")

    labs <- colheads
    if (!is.null(colnamesTexCmd))
      labs <- paste(sl, colnamesTexCmd, " ", labs, sep="")
                                        # DRW 12apr05.
    header <- NULL
    if(length(labs)) {
      if(!length(extracolheads)) {
        heads <- get2rowHeads(labs)
        colheads <- heads[[1]]
        if(any(heads[[2]] != ''))
          extracolheads <- heads[[2]]
      }
      
      if(multicol) ## SSJ 17nov03
        colheads <- paste(slmc1, cvbar, "}{", colheads, "}", sep="")
      
      header <- paste(colheads, collapse='&')
      if(length(extracolheads)) {
        extracolheads <- ifelse(extracolheads==''| extracolsize=='',
                                extracolheads,
                                paste('{',sl,extracolsize,' ',
                                      extracolheads,'}',sep=''))
        
        ## SSJ 17nov03 add | extracolsize=='' to avoid putting {\ } if you don't wont change size in second line title 
        if(multicol) ## SSJ 17nov03
          extracolheads <- ifelse(extracolheads=='',extracolheads,
                                  paste(slmc1,cvbar,'}{',extracolheads,'}',sep=''))
        else
          extracolheads <- ifelse(extracolheads=='',extracolheads,
                                  paste(extracolheads,sep=''))
        
        ##cat(eol," ", paste(c(if(length(rowname))'',extracolheads),collapse='&'),
        ##file=file, append=file!='') # 21jan03
        header <- paste(header, eol, paste(extracolheads, collapse='&'), sep='')
      }
    
      cat(header, eog, file=file, sep='', append=file!='') # 28apr03

      if(ctable)
        cat(midrule, file=file, append=file!='')
      else
        cat(midrule, file=file, append=file!='')
      ## eol was sl, sl  13dec02
    }
  }


  if(longtable) {
    if(missing(caption))
      cat(sl,"endhead\n",midrule,sl,"endfoot\n",sep="",
          file=file,append=file!='')
    else {
      cat(sl,"endfirsthead", sep="",file=file, append=file!='')
      cat(sl,"caption[]{\\em (continued)} ", eol,
          sep="",file=file, append=file!='')
      cat(midrule, sep="",file=file, append=file!='')
      cat(header, file=file, sep="&", append=file!='')
      cat(eog, midrule, sl, "endhead", '\n', midrule,
          sep="", file=file, append=file!='')
      if(length(insert.bottom)) {
        cat(paste(sl, 'multicolumn{', nc, '}{', "p{",sl,'linewidth}}{', 
                  insert.bottom, '}', eol, sep='', collapse='\n'),
                  sep="", file=file, append=file!='')
      }
    
      cat(sl,"endfoot\n", sep="",file=file, append=file!='')
      cat(sl,"label{", label, "}\n", sep="", file=file, append=file!='')
    }
  }

  { # individual lines, grouped if appropriate, longtable if appropriate
    if (length(n.rgroup)) {
      rg.end   <- cumsum(n.rgroup)
      rg.start <- rg.end-n.rgroup+1
      if(!length(rgroup)) {
        rgroup <- rep("",length(n.rgroup))
      } else {
        if (!is.null(rgroupTexCmd)) { # DRW 12apr05. This if block.
          rgroup <- paste("{",sl, rgroupTexCmd, " ", rgroup,"}",sep="") 
        } else {
          rgroup <- paste("{", rgroup,"}",sep="") 
        }
      }
      
      ##else rgroup <- paste("{",sl,"bf ",rgroup,"}",sep="") 
      seq.rgroup <- seq(along=n.rgroup)
    } else {
      seq.rgroup <- 1
      rg.end <- nr
      rg.start <- 1
    }

    linecnt <- 0
    for (j in seq.rgroup) {
      if (length(n.rgroup)) {
        if(longtable && linecnt>0 &&
           (linecnt+n.rgroup[j]+(n.rgroup[j]>1)) > lines.page) {
          cat(sl,"newpage\n", sep="",file=file, append=file!='')
          linecnt <- 0
        }
        
        cat(rgroup[j], rep("",nc-1), sep="&", file=file, append=file!='')
        cat(eol, sep="",file=file, append=file!='')
        ## eol was sl,sl 13dec02
        linecnt <- linecnt+1
      }

      ## Write the object (and it's formatting instructions)
      ## to the output.
      ## Loop through the rows of the object.
      for(i in rg.start[j]:rg.end[j]) {
        if (!length(n.rgroup)) {
          if(longtable && linecnt>0 && (linecnt+1 > lines.page)) {
            cat(sl,"newpage\n",sep="",file=file, append=file!='')
            linecnt <- 0						
          }
        }

        ## Loop through the columns of the object
        ## write each value (and it's format if there
        ## is one). 
        ## DRW 12apr05. This if/else block.
        if (!is.null(rcellTexCmds)) {
          num.cols <- ncol(cx)
          for (colNum in 1:num.cols) {
            cat(rcellTexCmds[i, colNum], " ", cx[i, colNum],
                file=file, append=file!='')
            if (colNum < num.cols)
              cat("&", file=file, append=file!='')
          }
        } else {
          ## Original code that writes object to output.
          cat(cx[i,], file=file, sep="&", append=file!='')
        }
        
        cat(if(i == rg.end[j] || (!ctable && !length(n.rgroup)))
              eog
            else if(i < rg.end[j])
              eol,
            sep="",file=file, append=file!='')
        
        ## eol was sl,sl  added if( ) 13dec02
        linecnt <- linecnt+1
      }  ## End of for loop that writes the object.

      if(length(n.rgroup) > j)
        cat(midrule, sep = "", file=file, append=file!='')
      else
        cat(bottomrule, sep="",file=file, append=file!='')
    }
  }

    if(do.end)
        cat(latex.end, file=file, sep="\n", append=file!='')
    
  sty <- c("longtable"[longtable], "here"[here], "dcolumn"[dcolumn],
           "ctable"[ctable], "booktabs"[booktabs],
           if(landscape && !ctable) "lscape")
  
  structure(list(file=file, style=sty), class='latex')
}
