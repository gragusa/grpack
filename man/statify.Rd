\name{statify}
\alias{statify}
\title{statify regression object}
\usage{
  statify(x,
    vcov. = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"),
    full = TRUE, mask, df = NULL, ...)
}
\arguments{
  \item{x}{a suitable regression object}

  \item{vcov.}{the variance to be used}

  \item{full}{if full output is needed}

  \item{mask}{coefficient to be exluded}

  \item{df}{degrees of freedom for the p-values}

  \item{...}{}
}
\value{
  Output
}
\description{
  "Statify" regression output
}
\details{
  Given a regression object out a stata-like table of
  results.
}
\author{
  Giuseppe M. Ragusa
}

