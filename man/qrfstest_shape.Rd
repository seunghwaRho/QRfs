\name{qrfstest_shape}
\alias{qrfstest_shape}
\title{Shape Tests of One Coefficient}
\description{
Tests that one coefficient is constant, or monotone, across quantiles (Chen et al. 2026, Section 4.4): the
homogeneity test, the monotonicity test over adjacent quantiles, or the monotonicity test over
every pair of quantiles.
}
\usage{
qrfstest_shape(arg_formula, arg_data, arg_variable, arg_har_method, arg_test,
               arg_alternative = c("decreasing", "increasing"),
               arg_tau_grid = NULL, arg_M = "test_optimal",
               arg_K = "test_optimal", arg_M_lambda = "mse_optimal_ar",
               arg_K_lambda = "mse_optimal_var", arg_alpha = 0.05,
               arg_kappa = NULL, arg_threshold_method = "theoretical",
               arg_var_fit = "gaussian",
               arg_var_stationary = "andrews_monahan", arg_coverage = 0.99,
               arg_hinv_method = "ker", arg_rq_method = "br",
               arg_seed = NULL)
}
\arguments{
  \item{arg_formula}{A formula, e.g. \code{y ~ x1 + x2}. Rows with missing values
      are dropped.}
  \item{arg_data}{A data frame holding the variables of the formula.}
  \item{arg_variable}{Coefficient name as it appears in the model matrix, e.g. \code{"NFCI"}.}
  \item{arg_har_method}{\code{"ker-har"} or \code{"os-har"}.}
  \item{arg_test}{\code{"homogeneity"} tests whether the coefficient is constant.
      It is two-sided over the adjacent differences. \code{"monotonicity_adjacent"}: one-sided over the adjacent
      differences, as in Chen et al. (2026, Section 5). \code{"monotonicity_allpair"}: one-sided over every pair of
      quantiles in the grid, following Chen et al. (2026, Section 4.4) with the sign convention below.}
  \item{arg_alternative}{Applicable to the monotonicity tests. \code{"decreasing"} \strong{(default)}
      tests the non-decreasing null. The \code{"increasing"} option tests the non-increasing null.
      Supplying this argument with the homogeneity test produces a warning, and the test proceeds.}
  \item{arg_tau_grid}{Quantile grid used for the test. When not specified (\code{NULL}),
      the default is \eqn{\{0.1,0.3,0.5,0.7,0.9\}} over \eqn{(0,1)}. Supply a vector such as
      \code{c(0.1, 0.5, 0.9)} to specify the grid directly. Alternatively, supply
      \code{list(range = c(0.1, 0.9), n = 5)} to generate a grid over a specified range.}
  \item{arg_alpha}{Significance level. Default \eqn{0.05}.}
  \item{arg_M}{Bandwidth \eqn{M} of the statistic under \code{"ker-har"}:
      \code{"test_optimal"} \strong{(default)}, \code{"mse_optimal_var"}, \code{"mse_optimal_ar"},
      \code{"simple"}, or a user-provided integer in \eqn{[1, T]}.}
  \item{arg_K}{Basis count \eqn{K} under \code{"os-har"}: \code{"test_optimal"} \strong{(default)},
      \code{"mse_optimal_var"}, \code{"mse_optimal_ar"}, \code{"simple"},
      or a user-provided integer in \eqn{[1,T-1]}.}
  \item{arg_M_lambda}{\eqn{M_\Lambda} used to estimate nuisance parameter \eqn{\Lambda} for kernel critical-value simulation:
      \code{"mse_optimal_ar"} \strong{(default)}, \code{"mse_optimal_var"},
      \code{"simple"}, or a user-provided integer.}
  \item{arg_K_lambda}{\eqn{K_\Lambda} used to estimate nuisance parameter \eqn{\Lambda} for OS critical-value simulation:
      \code{"mse_optimal_var"} \strong{(default)}, \code{"mse_optimal_ar"},
      \code{"simple"}, or a user-provided integer.}
  \item{arg_coverage}{Fraction of positive-eigenvalue variance retained in \eqn{\hat\Lambda}.
    Default \eqn{0.99}, permitted range \eqn{(0,0.99]}.}
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
For the coefficient \eqn{\beta_j(\tau)} and two grid quantiles \eqn{\tau_l<\tau_i}, define
\deqn{t_{il}=\frac{\sqrt{T}\,[\hat\beta_j(\tau_i)-\hat\beta_j(\tau_l)]}{\hat\sigma_j(\tau_i,\tau_l)},}
where \eqn{\hat\sigma_j} estimates the standard deviation of
\eqn{\sqrt{T}[\hat\beta_j(\tau_i)-\hat\beta_j(\tau_l)]}. The three tests differ in the pairs
compared and in the direction:
\itemize{
  \item \code{"homogeneity"}: \eqn{H_0}: \eqn{\beta_j(\tau)} is constant over \eqn{\mathcal T}.
    The two-sided statistic is \eqn{t^H=\max_{2\le i\le m}|t_{i,i-1}|} over the adjacent pairs.
  \item \code{"monotonicity_adjacent"}: \eqn{\max_{2\le i\le m}(-t_{i,i-1})} over the
    \eqn{m-1} adjacent grid pairs (Chen et al. 2026, Section 5).
  \item \code{"monotonicity_allpair"}: \eqn{\max_{1\le l<i\le m}(-t_{il})} over all
    \eqn{m(m-1)/2} grid pairs (Chen et al. 2026, Section 4.4).
}
Both monotonicity procedures test the same null over the quantile range:
\eqn{H_0:\beta_j(\tau_2)-\beta_j(\tau_1)\ge0} for all \eqn{\tau_1<\tau_2} in \eqn{\mathcal T}.
The formulas above use the default \code{"decreasing"} alternative.
\code{"increasing"} reverses the null inequality and uses \eqn{t} instead of \eqn{-t} in each maximum.
For \code{"decreasing"}, the library reports the negative of the minimum \eqn{t}-statistic
and rejects for large values.
}
\value{
A list of class \code{"qrfstest"} with the components

  \item{hypothesis}{The null hypothesis, as printed.}
  \item{coefficients}{The tested coefficient at each grid quantile.}
  \item{stat}{The maximum defined above.}
  \item{cv}{The critical value at level \code{alpha}.}
  \item{p_value}{The p-value.}
  \item{reject}{\eqn{1} if \code{stat > cv}, otherwise \eqn{0}.}
  \item{M_effective}{The bandwidth \eqn{M} of the statistic under \code{"ker-har"}, and \code{NA} under \code{"os-har"}.}
  \item{K_effective}{The basis count \eqn{K} of the statistic under \code{"os-har"}, and \code{NA} under \code{"ker-har"}.}
  \item{M_lambda, K_lambda}{The nuisance \eqn{M} or \eqn{K} used for the long-run covariance of the critical-value simulation. The one not used is \code{NA}.}
  \item{smoothing_rule}{\code{arg_M} or \code{arg_K} as given: the rule, or the number.}
  \item{tau_range}{The range of the null, as printed.}
  \item{q}{The number of pairs, \eqn{m-1} or \eqn{m(m-1)/2}.}
  \item{type}{\code{arg_test}.}
  \item{har_method}{\code{arg_har_method}.}
  \item{alpha}{The significance level.}
  \item{n_cv}{The number of simulated draws behind the critical value.}
  \item{nobs}{The number of observations \eqn{T}.}
  \item{tau}{The quantile grid.}
  \item{difference, t_ratio}{The tested differences and their \eqn{t}-statistics, both labelled by the quantile pair (\code{"0.3-0.1"}).}
}
\references{
Chen, K., Galvao, A. F., Rho, S., Vogelsang, T. J. and Yoon, J. (2026). Fixed-smoothing uniform inference for quantile regression. Working paper.
}
\seealso{
\code{\link{qrfstest_sup}}, \code{\link{qrfstest_stacked}}, \code{\link{print.qrfstest}}, \code{\link{summary.qrfstest}}
}
\examples{
\donttest{
data <- adrian
qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, data, "ggdp", "os-har",
               arg_test = "homogeneity", arg_seed = 1)

qrfstest_shape(
  ggdp_h1 ~ ggdp + NFCI, data, "ggdp", "os-har",
  arg_test = "homogeneity", arg_tau_grid = c(0.1, 0.25, 0.5, 0.75, 0.9),
  arg_seed = 1
)

qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, data, "NFCI", "ker-har",
               arg_test = "monotonicity_adjacent", arg_seed = 1)

qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, data, "NFCI", "ker-har",
               arg_test = "monotonicity_allpair", arg_alternative = "increasing",
               arg_seed = 1)
}
}
