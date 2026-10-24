\name{uniform_CI_band}
\alias{uniform_CI_band}
\title{Plot the Uniform Confidence Band}
\description{
Computes and plots the uniform confidence band with a single restriction.
}
\usage{
uniform_CI_band(arg_result, main = NULL, xlab = expression(tau), ylab = NULL,
          ylim = NULL, ...)
}
\arguments{
  \item{arg_result}{An object of class \code{"qrfstest"} from \code{qrfstest_sup} with one restriction.
      Any other object stops with a message.}
  \item{main}{The title. By default, it shows the confidence level and the HAR family.}
  \item{xlab}{The horizontal label. Default \eqn{\tau}.}
  \item{ylab}{The vertical label. When not specified (\code{NULL}), labels the tested coefficient combination.}
  \item{ylim}{The vertical limits. When not specified (\code{NULL}), includes the entire band and null value.}
  \item{\dots}{Passed on to \code{plot}.}
}
\value{
\code{band}, a data frame with \code{tau}, \code{estimate}, \code{lower} and
\code{upper}, invisibly. Called for the plot.
}
\seealso{
\code{\link{qrfstest_sup}}, \code{\link{plot.qrfstest}}
}
