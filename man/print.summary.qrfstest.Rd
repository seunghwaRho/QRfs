\name{print.summary.qrfstest}
\alias{print.summary.qrfstest}
\title{Print a Test Summary}
\description{
Prints a test summary.
}
\usage{
\method{print}{summary.qrfstest}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
  \item{x}{An object of class \code{"summary.qrfstest"} from \code{summary.qrfstest}.}
  \item{digits}{Significant digits of the printed numbers. Default
      \eqn{\max(3, \code{getOption("digits")} - 3)}.}
  \item{\dots}{Not used.}
}
\value{
Its argument, invisibly. Called for the printout.
}
\seealso{
\code{\link{summary.qrfstest}}
}
