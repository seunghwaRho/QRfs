\name{qrfs_example}
\alias{qrfs_example}
\title{Run the Worked Examples}
\description{
Runs and prints example tests using the Adrian, Boyarchenko and Giannone (2019) data.
}
\usage{
qrfs_example(arg_test = c("stacked", "sup", "shape"),
             arg_har_method = "ker-har", arg_var_fit = "gaussian",
             arg_seed = NULL)
}
\arguments{
  \item{arg_test}{Which example to run: \code{"stacked"}, \code{"sup"} or \code{"shape"}.
          The \code{"shape"} option runs the homogeneity test of \code{ggdp} and three
          monotonicity tests of \code{NFCI}.}
  \item{arg_har_method}{\code{"ker-har"} \strong{(default)} or \code{"os-har"}. Every command in
          the chosen example uses it.}
  \item{arg_var_fit}{\code{"gaussian"} \strong{(default)} or \code{"ols"}, passed to each
          test for its AR(1) or VAR(1) score approximation.}
  \item{arg_seed}{Seed passed to each test as its \code{arg_seed}. When not specified
          (\code{NULL}), the tests use the current random-number stream. Supply a seed to make
          the printed results reproducible. Each printed command can then be run on its own.}
}
\value{
Invisibly, a list of the \code{"qrfstest"} results, one per command run. The output is printed.
}
\seealso{
\code{\link{qrfstest_stacked}}, \code{\link{qrfstest_sup}}, \code{\link{qrfstest_shape}}, \code{\link{adrian}}
}
