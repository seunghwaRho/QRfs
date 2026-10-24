\name{adrian}
\alias{adrian}
\docType{data}
\title{Adrian, Boyarchenko and Giannone (2019) Data}
\description{
171 observations for one-quarter-ahead GDP growth regressions using the
Adrian, Boyarchenko and Giannone (2019) data. Forecast origins are 1973Q1--2015Q3.
}
\usage{data(adrian)}
\format{
GGDP is the annualized real GDP growth rate. The data frame contains \code{quarter}
(\code{"1973Q1"} and so on), \code{ggdp_h1} (the rate in quarter \eqn{t+1}),
\code{ggdp} (the rate in quarter \eqn{t}), and \code{NFCI}
(the Chicago Fed National Financial Conditions Index in quarter \eqn{t}).
}
\source{
Adrian, T., Boyarchenko, N. and Giannone, D. (2019). Vulnerable Growth. \emph{American Economic
Review} 109(4), 1263--1289, \doi{10.1257/aer.20161923}. Replication data \code{DataVulnerability.xls}, openICPSR 113169,
distributed under the Creative Commons Attribution 4.0 license (CC BY 4.0).
}
\seealso{
\code{\link{qrfs_example}}
}
\examples{
head(adrian, 3)
qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, adrian,
                 "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", "os-har")
}
