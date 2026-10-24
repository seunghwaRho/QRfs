\name{coef.qrfstest}
\alias{coef.qrfstest}
\title{Extract the Estimates}
\description{
Returns the estimates a result carries: the \eqn{p\times m} coefficient matrix of a stacked or sup test,
or the tested coefficient at each grid quantile for a shape test.
}
\usage{
\method{coef}{qrfstest}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{"qrfstest"}.}
  \item{\dots}{Ignored. Present for consistency with the \code{coef} generic.}
}
\value{
The \code{coefficients} field of the result.
}
\seealso{
\code{\link{qrfstest_stacked}}, \code{\link{qrfstest_sup}}, \code{\link{qrfstest_shape}}, \code{\link{qrhac}}
}
