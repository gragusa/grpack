\name{reg}
\alias{reg}
\title{lm}
\usage{
  reg(formula, data, subset, weights, na.action,
    method = "qr", model = TRUE, x = FALSE, y = FALSE,
    qr = TRUE, singular.ok = TRUE, contrasts = NULL,
    offset, cluster, ...)
}
\arguments{
  \item{formula}{an object of class \code{formula} (or one
  that can be coerced to that class): a symbolic
  description of the model to be fitted. See \code{lm} for
  details.}

  \item{data}{an optional data frame, list or environment
  (or object coercible by as.data.frame to a data frame)
  containing the variables in the model. If not found in
  data, the variables are taken from environment(formula),
  typically the environment from which lm is called.}

  \item{subset}{an optional vector specifying a subset of
  observations to be used in the fitting process.}

  \item{weights}{an optional vector of weights to be used
  in the fitting process. Should be NULL or a numeric
  vector. If non-NULL, weighted least squares is used with
  weights weights (that is, minimizing sum(w*e^2));
  otherwise ordinary least squares is used. See \code{lm}
  for "details".}

  \item{na.action}{a function which indicates what should
  happen when the data contain NAs. The default is set by
  the na.action setting of options, and is na.fail if that
  is unset. The "factory-fresh" default is na.omit. Another
  possible value is NULL, no action. Value na.exclude can
  be useful.}

  \item{method}{the method to be used; for fitting,
  currently only \code{method = "qr"} is supported;
  \code{method = "model.frame"} returns the model frame
  (the same as with \code{model = TRUE}, see below).}

  \item{model,}{y, x, qr logical.  logicals.  If ‘TRUE’ the
  corresponding components of the fit (the model frame, the
  model matrix, the response, the QR decomposition) are
  returned.}

  \item{singular.ok}{logical. If FALSE (the default in S
  but not in R) a singular fit is an error.}

  \item{contrasts}{an optional list. See the contrasts.arg
  of model.matrix.default.}

  \item{offset}{this can be used to specify an a priori
  known component to be included in the linear predictor
  during fitting. This should be NULL or a numeric vector
  of length equal to the number of cases. One or more
  offset terms can be included in the formula instead or as
  well, and if more than one are specified their sum is
  used. See \code{model.offset}}

  \item{cluster}{a factor specifying at which level
  clustering of standard error should take place.}

  \item{...}{additional arguments to be passed to the low
  level regression fitting functions (see \code{lm}).}
}
\value{
  A list similar to the one returned by \code{lm}.
}
\description{
  \code{reg} is used to fit linear models. It can be used
  to carry out regression, single stratum analysis of
  variance and analysis of covariance (although aov may
  provide a more convenient interface for these). It
  extends \code{lm} by allowing cluster standard error and
  by defining a summary method which uses by default
  heteroskedastic robust standard errors.
}
\details{
  Details are similar to \code{lm}.

  The only difference is that the factor at which level
  standard errors are clustered must be specified.
}
\author{
  Giuseppe Ragusa
}

