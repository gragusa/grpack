\name{wildboot}
\alias{wildboot}
\alias{wildboot.reg}
\title{wildboot}
\usage{
  wildboot(object, ...)

  \method{wildboot}{reg} (object, reps = 999, null,
    type = c("radamacher", "mtp", "mn1", "mn2"), ...)
}
\arguments{
  \item{a}{suitable regression object}

  \item{reps}{number of bootstrap replications}

  \item{null}{a named list containing the null hypothesis}

  \item{type}{type of weighting (see details)}

  \item{...}{other arguments}
}
\value{
  An object of class \code{wild.reg}

  \code{NULL}
}
\description{
  Wild Bootstrap for regression model.
}
\details{
  Calculate several wild bootstrapped t-statistics.
}
\examples{
data(CASchools)
CASchools <- transform(CASchools, testscore=(math+read)/2, str = students/teachers)


lm1 <- reg(testscore~str, cluster=county, data = CASchools)
wb1 <- wildboot(lm1, null = list(str=0))
}
\author{
  Giuseppe Ragusa
}

