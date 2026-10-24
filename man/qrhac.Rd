\name{qrhac}
\alias{qrhac}
\title{Wald Test with Small-b and Large-K Asymptotic Distribution}
\description{
Performs traditional heteroskedasticity and autocorrelation consistent (HAC) inference, i.e. small-\eqn{b}
and large-\eqn{K} based inference.
}
\usage{
qrhac(arg_formula, arg_data, arg_hypothesis, arg_har_method,
      arg_M = "mse_optimal_ar", arg_K = "mse_optimal_var",
      arg_alpha = 0.05, arg_kappa = NULL,
      arg_threshold_method = "theoretical",
      arg_var_fit = "gaussian",
      arg_var_stationary = "andrews_monahan",
      arg_hinv_method = "ker", arg_rq_method = "br")
}
\arguments{
  \item{arg_formula}{A formula, e.g. \code{y ~ x1 + x2}. Rows with missing values
      are dropped.}
  \item{arg_data}{A data frame holding the variables of the formula.}
  \item{arg_hypothesis}{One or more linear hypothesis strings, with case-sensitive coefficient
      names and the quantile in brackets, e.g. \code{"NFCI[0.05] = NFCI[0.5] = NFCI[0.95]"},
      as in \code{\link{qrfstest_stacked}}.}
  \item{arg_har_method}{\code{"ker-har"} or \code{"os-har"}.}
  \item{arg_alpha}{Significance level. Default \eqn{0.05}.}
  \item{arg_M}{Bandwidth \eqn{M} of the statistic under \code{"ker-har"}:
      \code{"mse_optimal_ar"} \strong{(default)}, \code{"mse_optimal_var"}, \code{"test_optimal"},
      \code{"simple"}, or a user-provided integer in \eqn{[1, T]}.}
  \item{arg_K}{Basis count \eqn{K} under \code{"os-har"}: \code{"mse_optimal_var"} \strong{(default)},
      \code{"mse_optimal_ar"}, \code{"test_optimal"}, \code{"simple"},
      or a user-provided integer, at least the number of restrictions and at most \eqn{T-1}.}
  \item{arg_kappa}{Size-inflation factor \eqn{\kappa>1} for test-optimal smoothing selection. When not specified
      (\code{NULL}), the default is \eqn{1.10} for \eqn{T<500} and \eqn{1.05} for \eqn{T\ge500}. Read only by
      \code{"test_optimal"}.}
  \item{arg_threshold_method}{Threshold for the Gaussian score approximation: \code{"theoretical"} \strong{(default)}
      uses \eqn{\Phi^{-1}(\tau)}. The \code{"empirical"} option uses \eqn{\hat q_\tau}, the empirical quantile of standardized residuals.}
  \item{arg_var_fit}{How the AR(1) or VAR(1) score approximation is estimated.
      \code{"gaussian"} \strong{(default)} combines the bivariate-normal approximation for
      indicator dependence with regressor dependence. \code{"ols"} fits the stacked scores
      by VAR(1) for the test-based and VAR-based rules, or each score column by AR(1) for the
      AR-based rule.
      \code{arg_threshold_method} is read only under \code{"gaussian"}.}
  \item{arg_var_stationary}{Stationarity bound for the AR(1) or VAR(1) approximation:
      \code{"andrews_monahan"} \strong{(\eqn{0.97}, default)}, \code{"sul_phillips_choi"}
      (\eqn{1-1/\sqrt{T}}), \code{"unrestricted"} (a fit that is then not
      stationary stops with a message), or a user-provided scalar in \eqn{(0,1)}.}
  \item{arg_hinv_method}{Hessian estimator from \code{quantreg}. The default \code{"ker"} uses Powell's
      kernel sandwich. The \code{"nid"} option uses the Hendricks--Koenker sparsity sandwich.}
  \item{arg_rq_method}{The solver of \code{quantreg::rq}. The default \code{"br"} uses the
      Barrodale--Roberts simplex. The \code{"fn"} option uses the Frisch--Newton interior point.
      For details please refer to \code{quantreg}.}
}
\details{
For \eqn{H_0:R\beta(\bar\tau_m)=r}, the Wald statistic \eqn{W} is computed as in
\code{\link{qrfstest_stacked}}. Under traditional, small-\eqn{b} or large-\eqn{K},
asymptotics \eqn{W\Rightarrow\chi^2_q}, where \eqn{q=\operatorname{rank}(R)}. The critical value is the
\eqn{(1-\alpha)} quantile of \eqn{\chi^2_q} and the p-value is the upper-tail probability at \eqn{W} under \eqn{\chi^2_q}.
The default rule for \eqn{M} or \eqn{K} of the statistic is MSE-optimal: \code{"mse_optimal_ar"} for \eqn{M}
and \code{"mse_optimal_var"} for \eqn{K}.
}
\value{
A list of class \code{"qrfstest"} with the components of
\code{\link{qrfstest_stacked}}, where now \code{cv_method} is
\code{"chisq"} and \code{n_cv} is \code{NA}.
}
\references{
Chen, K., Galvao, A. F., Rho, S., Vogelsang, T. J. and Yoon, J. (2026). Fixed-smoothing uniform inference for quantile regression. Working paper.

Galvao, A. F. and Yoon, J. (2024). HAC covariance matrix estimation in quantile regression. \emph{Journal of the American Statistical Association} 119, 2305--2316.
}
\seealso{
\code{\link{qrfstest_stacked}}, \code{\link{print.qrfstest}}, \code{\link{summary.qrfstest}}
}
\examples{
data <- adrian

qrhac(arg_formula = ggdp_h1 ~ ggdp + NFCI,
      arg_data = data,
      arg_hypothesis = "NFCI[0.05] = NFCI[0.5] = NFCI[0.95]",
      arg_har_method = "ker-har")

qrhac(ggdp_h1 ~ ggdp + NFCI, data,
      "NFCI[0.05] = NFCI[0.5] = NFCI[0.95]", "os-har")

# the test-based rule, or a bandwidth fixed by hand, in place of the MSE-optimal default
qrhac(ggdp_h1 ~ ggdp + NFCI, data,
      "NFCI[0.05] = NFCI[0.5] = NFCI[0.95]", "ker-har", arg_M = "test_optimal")

qrhac(ggdp_h1 ~ ggdp + NFCI, data,
      "NFCI[0.9] = 0", "ker-har")
}
