\name{summary.qrfstest}
\alias{summary.qrfstest}
\title{Summarize a Test Result}
\description{
Collects what a test result carries for printing in full: the estimates, the hypothesis, the test
statistic with its critical value and the limit that calibrates it (the fixed-\eqn{b} limit of
equation 18, the fixed-\eqn{K} limit of equation 19, the chi-square limit of \code{qrhac}, or Algorithm 1 of Chen et al. 2026), and the selected smoothing,
including the nuisance smoothing of a sup or shape test that \code{print.qrfstest} leaves out.
}
\usage{
\method{summary}{qrfstest}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{"qrfstest"}.}
  \item{\dots}{Not used.}
}
\value{
A list of class \code{"summary.qrfstest"} with one component, \code{result}, the test result. It is
printed by \code{\link{print.summary.qrfstest}}.
}
\references{
Chen, K., Galvao, A. F., Rho, S., Vogelsang, T. J. and Yoon, J. (2026). Fixed-smoothing uniform inference for quantile regression. Working paper.
}
\seealso{
\code{\link{print.summary.qrfstest}}, \code{\link{print.qrfstest}}
}
