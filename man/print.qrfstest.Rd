\name{print.qrfstest}
\alias{print.qrfstest}
\title{Print a Test Result}
\description{
Prints the null hypothesis, test statistic, critical value, p-value, sample size and selected
bandwidth fraction \eqn{b} or basis count \eqn{K}.
}
\usage{
\method{print}{qrfstest}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{"qrfstest"}, from any of the test functions.}
  \item{\dots}{Ignored. Present for consistency with the \code{print} generic.}
}
\details{
All tests print \code{H0}. Sup and shape tests state the null over the quantile range,
with the numerical grid on the next line.
Sup and shape tests also print their individual statistics. A sup test with a single
restriction prints them as a table with the uniform confidence band. The selected \eqn{b} or \eqn{K}
is labelled \code{data dependent} or \code{user provided}. A \code{qrhac} result is headed by its
method, small-\eqn{b} Bartlett kernel HAC or large-\eqn{K} orthonormal series HAR, in place of \code{HAR}.
}
\value{
Its argument, invisibly.
}
\seealso{
\code{\link{summary.qrfstest}}, \code{\link{qrfstest_stacked}}, \code{\link{qrfstest_sup}}, \code{\link{qrfstest_shape}}, \code{\link{qrhac}}
}
