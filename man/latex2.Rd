\name{latex2}
\alias{latex2}
\title{Latex tables.}
\usage{
  latex2(object,
    title = first.word(deparse(substitute(object))),
    file = paste(title, ".tex", sep = ""), append = FALSE,
    do.begin = TRUE, do.end = TRUE, label = title,
    rowlabel = title, rowlabel.just = "l", cgroup = NULL,
    n.cgroup = NULL, rgroup = NULL, n.rgroup = NULL,
    cgroupTexCmd = "bfseries", rgroupTexCmd = "bfseries",
    rownamesTexCmd = NULL, colnamesTexCmd = NULL,
    cellTexCmds = NULL, rowname,
    cgroup.just = rep("c", length(n.cgroup)),
    colheads = NULL, extracolheads = NULL,
    extracolsize = "scriptsize", dcolumn = FALSE,
    tight = FALSE, numeric.dollar = !dcolumn, cdot = FALSE,
    longtable = FALSE, draft.longtable = TRUE,
    ctable = FALSE, booktabs = FALSE, table.env = TRUE,
    here = FALSE, lines.page = 40, caption = NULL,
    caption.lot = NULL, caption.loc = c("top", "bottom"),
    double.slash = FALSE, vbar = FALSE,
    collabel.just = rep("c", nc), na.blank = TRUE,
    insert.bottom = NULL,
    first.hline.double = !(booktabs | ctable),
    where = "!tbp", size = NULL,
    center = c("center", "centering", "none"),
    landscape = FALSE, multicol = TRUE,
    math.row.names = FALSE, math.col.names = FALSE, ...)
}
\arguments{
  \item{do.begin}{Should the latex code at the beginning of
  table be printed?}

  \item{do.end}{Should the latex code at the end of table
  be printed?}

  \item{tight}{when using dcolumn tighten columns width.}
}
\value{
  See \code{\link{latex}} of package \cod{Hmisc}.
}
\description{
  Exactly like latex of \code{Hmisc}, but fix a bug when
  using \code{ctable=TRUE} and add three arguments:
  \cod{do.begin}, \code{do.end}, \code{tight}.
}
\details{
  .. content for \details{} ..
}
\author{
  Giuseppe Ragusa
}
