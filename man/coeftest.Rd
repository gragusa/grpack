\name{coeftest.reg}
\alias{coeftest.reg}
\title{Testing Estimated Coefficients}
\usage{
  coeftest.reg(x,
    vcov. = c("HC1", "const", "HC", "HC0", "HC2", "HC3", "HC4", "HC4m", "HC5", "HAC"),
    df = NULL)
}
\arguments{
  \item{x}{a \code{reg} object}

  \item{vcov}{a covariance type}

  \item{df}{the degrees of freedom to be used. If this is a
  finite positive number a t test with df degrees of
  freedom is performed. In all other cases, a z test (using
  a normal approximation) is performed. If the \code{reg}
  object has a \code{cluster} component and df is
  \code{NULL} a t test with G-1 degrees of freedom is
  performed.}

  \item{...}{other arguments}
}
\description{
  \code{coeftest} is a method for performing z and
  (quasi-)t tests of estimated coefficients through
  \code{reg}
}
\details{
  Details
}

