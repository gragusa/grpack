\name{bfgs}
\alias{bfgs}
\title{bfgs}
\usage{
  bfgs(x0, f_, g_, min.x = rep(-Inf, length(x0)),
    max.x = rep(Inf, length(x0)), prec = 1e-05,
    verbose = FALSE)
}
\arguments{
  \item{x0}{Tthe initial solution guess}

  \item{f_}{The function to be minimized}

  \item{g_}{The gradient of f_}

  \item{min.x}{lower bounds}

  \item{max.x}{upper bounds}

  \item{prec}{}

  \item{verbose}{if TRUE diplay iteration nformation}
}
\value{
  A list

  Apart from returning the maximizer x, it also returns the
  function value f(x) the gradient f'(x), and the inverse
  hessian f''(x)^{-1}.
}
\description{
  Implements the Broyden-Fletcher-Goldfarb-Shanno algorithm
  for function minimization.
}
\details{
  The function attempts to find

  argmin_{x s.t. min.x < x < max.x coord-wise} f(x)
}
\author{
  Giuseppe Ragusa
}

