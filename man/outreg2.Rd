\name{outreg2}
\alias{outreg2}
\title{Arrange regression outputs into an illustrative table}
\usage{
  outreg2(..., model.list, model.names = NULL,
    order = c("lr", "rl", "longest"), omitcoef = NULL,
    omit.model = NULL, coef.names = NULL,
    cons.below = TRUE, stars = "stata", lev = 0.05,
    additional.rows = NULL,
    title = first.word(deparse(substitute(object))),
    file = "", append = FALSE, label = title,
    rowlabel = "", rowlabel.just = "l", cgroup = NULL,
    n.cgroup = NULL, rgroup = NULL, n.rgroup = NULL,
    cgroupTexCmd = NULL, rgroupTexCmd = NULL,
    rownamesTexCmd = NULL, colnamesTexCmd = NULL,
    cellTexCmds = NULL, rowname,
    cgroup.just = rep("c", sum(n.cgroup)), colheads = NULL,
    extracolheads = NULL, extracolsize = "scriptsize",
    dcolumn = TRUE, tight = TRUE,
    numeric.dollar = !dcolumn, cdot = FALSE,
    longtable = FALSE, draft.longtable = TRUE,
    ctable = FALSE, booktabs = FALSE, table.env = FALSE,
    here = FALSE, lines.page = 40, caption = NULL,
    caption.lot = NULL, caption.loc = c("top", "bottom"),
    double.slash = FALSE, vbar = FALSE,
    collabel.just = rep("c", nc), na.blank = TRUE,
    insert.bottom = NULL, do.begin = TRUE, do.end = TRUE,
    first.hline.double = !(booktabs | ctable),
    where = "!tbp", size = NULL,
    center = c("center", "centering", "none"),
    landscape = FALSE, multicol = TRUE,
    math.row.names = TRUE, math.col.names = FALSE,
    rowcolors = NULL)
}
\arguments{
  \item{...}{the statistical models.}

  \item{model.list}{the statistical models as a list. If
  model list is missing, then the models are assumed to be
  passed by \code{...}.}

  \item{model.names}{Optional vector of names to use as
  column headings in the table. If more models than names
  are supplied, unnamed models are numbered (starting at
  one more than the number of names).}

  \item{order}{Determines the order in which terms (rows)
  are included in the output when more than one model
  (column) is present. “lr” and “rl” take the order of
  terms from the first or last (leftmost or rightmost)
  model and appends new terms as they are encountered.
  “longest” uses the order of terms in the model with the
  most \code{terms.Default = “lr”}.}

  \item{omitcoef}{An optional integer or character vector
  of coefficient indices, or an expression involving
  coefnames that evaluates to integer or character, of rows
  to exclude from the output. See details.}

  \item{omit.model}{An optional integer or character vector
  of model indices (in the order they are entred) to be
  excluded from the output. Useful if
  \code{order="largest"}, and the largest model is not to
  be shown. See details.}

  \item{coef.names}{An optional vector of names for
  coefficients. It is recommended to establish the
  \code{omitcoef} and \code{order} settings with automatic
  symbolic naming before supplying a vector of “pretty”
  variable names. If automatic symbolic naming is used,
  names are taken from the variables in the models and
  “sanitized” for latex. If \code{coef.names} are supplied,
  they must be valid latex, with double-backslash escape
  characters.}

  \item{cons.below}{}

  \item{stars}{how statistical significance “stars”, either
  "stata", "default", or 1. "stata" is based on the stata
  default and gives three stars.  "default" uses the
  \code{R} default, not to be confused with the function's
  (perhaps confusing) Default="stata"}

  \item{lev}{When \code{stars=1}, what level should be used
  for the test to reject statistical insignificance and
  bestow the glittering star? Disable decoration entirely
  by specifying \code{lev=0}. \code{Default=.05}.}

  \item{additional.rows}{A matrix with number of column
  equal to the number of included models for additional
  annotation of the table. Useful for fixed effects,
  country effects, etc. specifications.}

  \item{title}{}

  \item{file}{name of the file to create. The default is
  not to create a file and to have the generated LaTeX code
  just printed to standard output. This is especially
  useful when running under Sweave in R using its
  results=tex tag, to save having to manage many small
  external files. When file="", latex keeps track of LaTeX
  styles that are called for by creating or modifying an
  object \code{latexStyles} (in \code{.GlobalTemp} in R or
  in frame 0 in S-Plus). \code{latexStyles} is a vector
  containing the base names of all the unique LaTeX styles
  called for so far in the current session. See the end of
  the examples section for a way to use this object to good
  effect. For dvips, file is the name of an output
  postscript file.}

  \item{append}{defaults to FALSE. Set to TRUE to append
  output to an existing file.}

  \item{label}{a text string representing a symbolic label
  for the table for referencing in the LaTeX \code{\label}
  and \code{\ref} commands. \code{label} is only used if
  caption is given.}

  \item{rowlabel}{}

  \item{rowlabel.just}{}

  \item{cgroup}{}

  \item{n.cgroup}{}

  \item{rgroup}{}

  \item{n.rgroup}{}

  \item{cgroupTexCmd}{}

  \item{rgroupTexCmd}{}

  \item{rownamesTexCmd}{}

  \item{colnamesTexCmd}{}

  \item{cellTexCmds}{}

  \item{rowname}{}

  \item{cgroup.just}{}

  \item{colheads}{}

  \item{extracolheads}{}

  \item{extracolsize}{}

  \item{dcolumn}{}

  \item{tight}{}

  \item{numeric.dollar}{}

  \item{cdot}{}

  \item{longtable}{}

  \item{draft.longtable}{}

  \item{ctable}{}

  \item{booktabs}{}

  \item{table.env}{}

  \item{here}{}

  \item{lines.page}{}

  \item{caption}{}

  \item{caption.lot}{}

  \item{caption.loc}{}

  \item{double.slash}{}

  \item{vbar}{}

  \item{collabel.just}{}

  \item{na.blank}{}

  \item{insert.bottom}{}

  \item{do.begin}{Should the latex code at the beginning of
  table be printed? (Default: TRUE)}

  \item{do.end}{Should the latex code at the end of table
  be printed? (Default: TRUE)}

  \item{first.hline.double}{}

  \item{where}{}

  \item{size}{}

  \item{center}{}

  \item{landscape}{}

  \item{multicol}{}

  \item{math.row.names}{}

  \item{math.col.names}{}

  \item{rowcolors}{This command is useful for inserting
  color command. For instance, if
  \code{rowcolors=rowcolors{1}{green}{pink}}, then the rows
  are subsequently coloured in green and in pink.}
}
\value{
  A character vector containing the LaTeX table.
}
\description{
  Function to arrange regression outputs into an
  well-formatted LaTeX illustrative tables of multiple
  statistical models.
}
\details{
  \code{outreg2} borrows heavily from the \code{apsrtable}
  package. It preserves its flexibility in that many
  statistical models can be easily considered. It extends
  it by allowing finer control on the LaTeX output. This
  flexibility is obtained by using a modified version of
  the finction \code{latex} of the \code{Hmisc} package.
}
\author{
  Giuseppe Ragusa
}

