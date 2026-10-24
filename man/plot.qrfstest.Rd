\name{plot.qrfstest}
\alias{plot.qrfstest}
\title{Plot a Test Result}
\description{
Draws the uniform confidence band of a sup test with a single restriction by calling
\code{\link{uniform_CI_band}}.
}
\usage{
\method{plot}{qrfstest}(x, y, ...)
}
\arguments{
  \item{x}{An object of class \code{"qrfstest"} from \code{qrfstest_sup} with one restriction.
      Any other object stops with a message.}
  \item{y}{Not used. The generic is \code{plot(x, y, ...)}.}
  \item{\dots}{Passed on to \code{\link{uniform_CI_band}}.}
}
\value{
Its argument, invisibly. Called for the plot.
}
\seealso{
\code{\link{uniform_CI_band}}, \code{\link{qrfstest_sup}}
}
