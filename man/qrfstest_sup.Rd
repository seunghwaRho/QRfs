\name{qrfstest_sup}
\alias{qrfstest_sup}
\title{Sup Test over a Quantile Range}
\description{
Performs a uniform-in-\eqn{\tau} hypothesis test over a range of quantiles, with a simulated critical value and p-value (Chen et al. 2026, Section 4.1).
}
\usage{
qrfstest_sup(arg_formula, arg_data, arg_hypothesis, arg_har_method,
             arg_tau_grid = NULL, arg_M = "test_optimal",
             arg_K = "test_optimal", arg_M_lambda = "mse_optimal_ar",
             arg_K_lambda = "mse_optimal_var", arg_alpha = 0.05,
             arg_kappa = NULL, arg_threshold_method = "theoretical",
             arg_var_fit = "gaussian", arg_var_stationary = "andrews_monahan",
             arg_coverage = 0.99, arg_hinv_method = "ker",
             arg_rq_method = "br", arg_band_plot = FALSE, arg_seed = NULL)
}
\arguments{
  \item{arg_formula}{A formula, e.g. \code{y ~ x1 + x2}. Rows with missing values
      are dropped.}
  \item{arg_data}{A data frame holding the variables of the formula.}
  \item{arg_hypothesis}{Linear restrictions without quantile brackets, e.g. \code{"NFCI = 0"} or
      \code{c("NFCI = 0", "ggdp = 0")}.}
  \item{arg_har_method}{\code{"ker-har"} or \code{"os-har"}.}
  \item{arg_tau_grid}{Quantile grid used to discretize the range. When not specified (\code{NULL}),
      the default is \eqn{\{0.1,0.3,0.5,0.7,0.9\}} over \eqn{(0,1)}. Supply a vector such as
      \code{c(0.1, 0.5, 0.9)} to specify the grid directly. Alternatively, supply
      \code{list(range = c(0.1, 0.9), n = 5)} to generate a grid over a specified range.}
  \item{arg_band_plot}{\code{FALSE} \strong{(default)} or \code{TRUE}: draw the uniform confidence
      band when the call returns. Only a single restriction has a band. Otherwise a warning is given
      and nothing is drawn. The band can also be drawn later with \code{\link{uniform_CI_band}} or \code{plot()} on the result.}
  \item{arg_alpha}{Significance level. Default \eqn{0.05}.}
  \item{arg_M}{Bandwidth \eqn{M} of the statistic under \code{"ker-har"}:
      \code{"test_optimal"} \strong{(default)}, \code{"mse_optimal_var"}, \code{"mse_optimal_ar"},
      \code{"simple"}, or a user-provided integer in \eqn{[1, T]}.}
  \item{arg_K}{Basis count \eqn{K} under \code{"os-har"}: \code{"test_optimal"} \strong{(default)},
      \code{"mse_optimal_var"}, \code{"mse_optimal_ar"}, \code{"simple"},
      or a user-provided integer, at least the number of restrictions at each quantile and at most \eqn{T-1}.}
  \item{arg_M_lambda}{\eqn{M_\Lambda} used to estimate nuisance parameter \eqn{\Lambda} for kernel critical-value simulation:
      \code{"mse_optimal_ar"} \strong{(default)}, \code{"mse_optimal_var"},
      \code{"simple"}, or a user-provided integer.}
  \item{arg_K_lambda}{\eqn{K_\Lambda} used to estimate nuisance parameter \eqn{\Lambda} for OS critical-value simulation:
      \code{"mse_optimal_var"} \strong{(default)}, \code{"mse_optimal_ar"},
      \code{"simple"}, or a user-provided integer.}
  \item{arg_coverage}{The fraction of the positive-eigenvalue variance of \eqn{\hat\Lambda} kept by the
      Karhunen--Lo\`eve truncation of Chen et al. (2026, Algorithm 1). Default \eqn{0.99}, permitted range \eqn{(0, 0.99]}.}
  \item{arg_kappa}{Size-inflation factor \eqn{\kappa>1} for test-optimal smoothing selection. When not specified
      (\code{NULL}), the default is \eqn{1.10} for \eqn{T<500} and \eqn{1.05} for \eqn{T\ge500}.}
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
  \item{arg_seed}{Seed passed to \code{set.seed()} before the test. When not specified
      (\code{NULL}), the current random-number stream is used.}
}
\details{
Tests the hypotheses in Chen et al. (2026, Section 4.1):
\deqn{H_0: r(\beta(\tau)) = 0 \mbox{ for all } \tau \in \mathcal{T}, \qquad H_A: r(\beta(\tau)) \ne 0 \mbox{ for some } \tau \in \mathcal{T}.}
The library accepts linear restrictions, \eqn{r(\beta(\tau))=R_0\beta(\tau)-r_0},
and uses a quantile grid to approximate \eqn{\mathcal T}.
At each quantile, let \eqn{d_i=R_0\hat\beta(\tau_i)-r_0}. The statistic is
\deqn{W_{\mathrm{sup}}=\max_{1\le i\le m}W(\tau_i),\qquad
 W(\tau_i)=T d_i'[R_0\hat\Sigma(\tau_i)R_0']^{-1}d_i,}
where \eqn{\hat\Sigma(\tau_i)=\hat D(\tau_i)^{-1}\hat\Lambda(\tau_i,\tau_i)\hat D(\tau_i)^{-1}}.
Reject for large values of \eqn{W_{\mathrm{sup}}}. Critical values are simulated using
the estimated long-run covariance across quantiles (Chen et al. 2026, Algorithm 1).
The default rule for \eqn{M} or \eqn{K} is \code{"test_optimal"}.

For a single restriction (\eqn{q_0=1}) the result also carries the uniform confidence band of
Chen et al. (2026, Section 4.3, Algorithm 2). The sup-Wald critical value is the square of the sup-\eqn{|t|} one, so
\deqn{I_\alpha(\mathcal T)=\bigl[R_0\hat\beta(\tau)-\sqrt{\mathrm{cv}}\,\hat\sigma(\tau),\;
                          R_0\hat\beta(\tau)+\sqrt{\mathrm{cv}}\,\hat\sigma(\tau)\bigr],\qquad \tau\in\mathcal T,}
where \eqn{\hat\sigma(\tau)^2=R_0\hat\Sigma(\tau)R_0'/T} is the variance of the tested combination at
the test's own \eqn{M} or \eqn{K}.
}
\value{
A list of class \code{"qrfstest"} with the components

  \item{hypothesis}{The null hypothesis, as printed.}
  \item{coefficients}{The \eqn{p\times m} matrix of quantile regression estimates, one column per grid quantile.}
  \item{stat}{The sup-Wald statistic, the maximum of \code{wald_tau}.}
  \item{cv}{The critical value at level \code{alpha}.}
  \item{p_value}{The p-value.}
  \item{reject}{\eqn{1} if \code{stat > cv}, otherwise \eqn{0}.}
  \item{M_effective}{The bandwidth \eqn{M} of the statistic under \code{"ker-har"}, and \code{NA} under \code{"os-har"}.}
  \item{K_effective}{The basis count \eqn{K} of the statistic under \code{"os-har"}, and \code{NA} under \code{"ker-har"}.}
  \item{M_lambda, K_lambda}{The nuisance \eqn{M} or \eqn{K} used for the long-run covariance of the critical-value simulation. The one not used is \code{NA}.}
  \item{q}{The number of restrictions at each quantile.}
  \item{type}{\code{"sup"}.}
  \item{har_method}{\code{arg_har_method}.}
  \item{cv_method}{\code{"nonpivotal-sim"}.}
  \item{alpha}{The significance level.}
  \item{n_cv}{The number of simulated draws behind the critical value.}
  \item{nobs}{The number of observations \eqn{T}.}
  \item{tau}{The quantile grid.}
  \item{smoothing_rule}{\code{arg_M} or \code{arg_K} as given: the rule, or the number.}
  \item{tau_range}{The range of the null, as printed.}
  \item{R, null_val}{The restriction matrix and the null values of the hypothesis.}
  \item{wald_tau}{The \eqn{m} per-quantile Wald statistics.}
  \item{band}{For a single restriction, a data frame with \code{tau}, \code{estimate}, \code{lower} and \code{upper}, the uniform confidence band, and otherwise \code{NULL}.}
}
\references{
Chen, K., Galvao, A. F., Rho, S., Vogelsang, T. J. and Yoon, J. (2026). Fixed-smoothing uniform inference for quantile regression. Working paper.
}
\seealso{
\code{\link{uniform_CI_band}}, \code{\link{plot.qrfstest}}, \code{\link{qrfstest_stacked}}, \code{\link{qrfstest_shape}}, \code{\link{print.qrfstest}}, \code{\link{summary.qrfstest}}
}
\examples{
\donttest{
# nothing specified: the default grid discretizes the whole range
data <- adrian
qrfstest_sup(ggdp_h1 ~ ggdp + NFCI, data, "NFCI = 0", "ker-har",
             arg_seed = 1)

# the same null confined to the lower quantiles
qrfstest_sup(ggdp_h1 ~ ggdp + NFCI, data, "NFCI = 0", "ker-har",
             arg_tau_grid = c(0.05, 0.10, 0.15, 0.20, 0.25), arg_seed = 1)

# With the series. The band and the terms show the rejection comes from tau = 0.9 alone.
qrfstest_sup(ggdp_h1 ~ ggdp + NFCI, data, "ggdp = 0", "os-har",
             arg_seed = 1)
}
}
