\name{qrfstest_stacked}
\alias{qrfstest_stacked}
\title{Stack-Wald Test at Specified Quantiles}
\description{
Computes the stack-Wald statistic for joint linear restrictions at specified quantiles, its critical value and p-value (Chen et al. 2026, Section 4.2).
}
\usage{
qrfstest_stacked(arg_formula, arg_data, arg_hypothesis, arg_har_method,
                 arg_M = "test_optimal", arg_K = "test_optimal",
                 arg_alpha = 0.05, arg_kappa = NULL,
                 arg_threshold_method = "theoretical",
                 arg_var_fit = "gaussian",
                 arg_var_stationary = "andrews_monahan",
                 arg_hinv_method = "ker", arg_rq_method = "br",
                 arg_seed = NULL)
}
\arguments{
  \item{arg_formula}{A formula, e.g. \code{y ~ x1 + x2}. Rows with missing values
      are dropped.}
  \item{arg_data}{A data frame holding the variables of the formula.}
  \item{arg_hypothesis}{One or more linear hypothesis strings, with case-sensitive coefficient
      names and the quantile in brackets, e.g. \code{"NFCI[0.1] = NFCI[0.5]"},
      \code{"NFCI[0.1] = NFCI[0.5] = NFCI[0.9]"} or
      \code{"2*NFCI[0.25] + NFCI[0.75] = 1"}.
      To pass two strings together, use
      \code{c("NFCI[0.1] = NFCI[0.5] = NFCI[0.9]",}
      \code{"ggdp[0.1] = ggdp[0.9]")}.}
  \item{arg_har_method}{\code{"ker-har"} or \code{"os-har"}.}
  \item{arg_alpha}{Significance level. Default \eqn{0.05}.}
  \item{arg_M}{Bandwidth \eqn{M} of the statistic under \code{"ker-har"}:
      \code{"test_optimal"} \strong{(default)}, \code{"mse_optimal_var"}, \code{"mse_optimal_ar"},
      \code{"simple"}, or a user-provided integer in \eqn{[1, T]}.}
  \item{arg_K}{Basis count \eqn{K} under \code{"os-har"}: \code{"test_optimal"} \strong{(default)},
      \code{"mse_optimal_var"}, \code{"mse_optimal_ar"}, \code{"simple"},
      or a user-provided integer, at least the number of restrictions and at most \eqn{T-1}.}
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
For \eqn{H_0:R\beta(\bar\tau_m)=r}, the statistic is
\deqn{W=T(R\hat\beta-r)'[R\hat D^{-1}\hat\Lambda\hat D^{-1}R']^{-1}(R\hat\beta-r),}
where \eqn{\hat\beta} stacks the coefficient estimates and \eqn{\hat D} is block diagonal
across quantiles (Chen et al. 2026, equation 16). Reject for large values of \eqn{W}.

The Bartlett critical value is simulated from the pivotal fixed-\eqn{b} limit in Chen et al. (2026, equation 18).
The OS critical value follows from
\eqn{W^{\mathrm{OS}}\Rightarrow\frac{qK}{K-q+1}F_{q,K-q+1}}, where \eqn{q=\operatorname{rank}(R)}
(Chen et al. 2026, equation 19). The default rule for \eqn{M} or \eqn{K} is \code{"test_optimal"}.
}
\value{
A list of class \code{"qrfstest"} with the components

  \item{hypothesis}{The null hypothesis, as printed.}
  \item{coefficients}{The \eqn{p\times m} matrix of quantile regression estimates, one column per quantile.}
  \item{stat}{The stack-Wald statistic \eqn{W}.}
  \item{cv}{The critical value at level \code{alpha}.}
  \item{p_value}{The p-value.}
  \item{reject}{\eqn{1} if \code{stat > cv}, otherwise \eqn{0}.}
  \item{M_effective}{The bandwidth \eqn{M} of the statistic under \code{"ker-har"}, and \code{NA} under \code{"os-har"}.}
  \item{K_effective}{The basis count \eqn{K} of the statistic under \code{"os-har"}, and \code{NA} under \code{"ker-har"}.}
  \item{M_lambda, K_lambda}{\code{NA}: the stacked test estimates no nuisance covariance.}
  \item{q}{The number of restrictions.}
  \item{type}{\code{"stacked"}.}
  \item{har_method}{\code{arg_har_method}.}
  \item{cv_method}{\code{"pivotal"}.}
  \item{alpha}{The significance level.}
  \item{n_cv}{The number of simulated draws behind the critical value under \code{"ker-har"}, and \code{NA} under \code{"os-har"}, where the critical value is in closed form.}
  \item{nobs}{The number of observations \eqn{T}.}
  \item{tau}{The quantiles named in the hypothesis, sorted.}
  \item{smoothing_rule}{\code{arg_M} or \code{arg_K} as given: the rule, or the number.}
  \item{tau_range}{\code{NULL}.}
  \item{R, null_val}{The restriction matrix and the null values of the hypothesis.}
}
\references{
Chen, K., Galvao, A. F., Rho, S., Vogelsang, T. J. and Yoon, J. (2026). Fixed-smoothing uniform inference for quantile regression. Working paper.
}
\seealso{
\code{\link{qrhac}}, \code{\link{qrfstest_sup}}, \code{\link{qrfstest_shape}}, \code{\link{print.qrfstest}}, \code{\link{summary.qrfstest}}, \code{\link{coef.qrfstest}}
}
\examples{
data <- adrian

qrfstest_stacked(arg_formula = ggdp_h1 ~ ggdp + NFCI,
                 arg_data = data,
                 arg_hypothesis = "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]",
                 arg_har_method = "ker-har",
                 arg_seed = 1)

qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, data,
                 "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", "os-har")

# Two hypotheses in one test. The grid is the union of the quantiles they specify.
qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, data,
  c("NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", "ggdp[0.1] = ggdp[0.9]"), "os-har")

# smoothing fixed by hand, and a different level
qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, data,
                 "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", "os-har",
                 arg_K = 40, arg_alpha = 0.10)

# a larger kappa allows more size distortion, so the rule takes a larger K
for (k in list(NULL, 1.05, 1.50)) {
  r <- qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, data,
                        "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", "os-har",
                        arg_kappa = k)
  cat(sprintf("arg_kappa = \%-4s  K = \%3d   W = \%.4f   cv = \%.4f\n",
              if (is.null(k)) "NULL" else format(k), r$K_effective, r$stat, r$cv))
}
}
