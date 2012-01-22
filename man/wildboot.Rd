\name{wildboot}
\alias{wildboot}
\alias{wildboot.reg}
\title{wildboot}
\usage{
  wildboot(x, ...)

  \method{wildboot}{reg} (obj, reps = 999, null,
    type = c("radamacher", "mtp", "mn1", "mn2"))
}
\arguments{
  \item{x}{an object}

  \item{...}{other arguments}
}
\value{
  A list with several components

  \code{NULL}
}
\description{
  Wild Bootstrap for regression model.
}
\details{
  Calculate several wild bootstrapped quantities.
}
\author{
  Giuseppe Ragusa
}

