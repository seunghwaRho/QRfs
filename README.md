# QRfs

QRfs implements fixed-smoothing (fixed-*b* and fixed-*K*) methods for heteroskedasticity and
autocorrelation robust (HAR) inference in time-series quantile regression, following Chen, Galvao,
Rho, Vogelsang and Yoon (2026), *Fixed-smoothing Uniform Inference for Quantile Regression*.
Accounting for dependence both over time and across quantiles, it provides stack-Wald tests of
joint linear restrictions at specified quantiles, uniform sup tests and confidence bands over a
quantile range, and tests of coefficient homogeneity and monotonicity.

## Installation

```r
# install.packages("remotes")
remotes::install_github("seunghwaRho/QRfs")
```

## Examples

```r
library(QRfs)
```

The examples use the included `adrian` data from Adrian, Boyarchenko and Giannone (2019),
regressing one-quarter-ahead GDP growth (`ggdp_h1`) on current GDP growth (`ggdp`) and the
Chicago Fed National Financial Conditions Index (`NFCI`).

### Uniform sup test over a quantile range

```r
sup <- qrfstest_sup(ggdp_h1 ~ ggdp + NFCI, adrian,
                    "NFCI = 0", "ker-har",
                    arg_tau_grid = list(range = c(0, 0.25), n = 5))
uniform_CI_band(sup)
```

`qrfstest_sup()` performs a uniform-in-τ hypothesis test over a range of quantiles $\mathcal{T}$
with the sup-Wald statistic. The example above tests $H_0: \beta_{NFCI}(\tau) = 0$ for all
$\tau \in (0, 0.25]$, discretized at five grid points. `uniform_CI_band()` provides the uniform
confidence band.
`arg_har_method = "ker-har"` chooses fixed-*b* inference with the Bartlett kernel, `"os-har"`
fixed-*K* inference with the DCT-II orthonormal series.

### Stack-Wald test at specified quantiles

```r
qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, adrian,
                 "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", "os-har")
```

`qrfstest_stacked()` performs a test of joint linear restrictions at specified quantiles. The
example above tests $H_0: \beta_{NFCI}(0.1) = \beta_{NFCI}(0.5) = \beta_{NFCI}(0.9)$ with the
stack-Wald statistic.

### Homogeneity test

```r
qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, adrian, "NFCI", "ker-har",
               arg_test = "homogeneity")
qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, adrian, "NFCI", "ker-har",
               arg_test = "homogeneity",
               arg_tau_grid = c(0.1, 0.25, 0.5, 0.75, 0.9))
```

The homogeneity test assesses whether the NFCI coefficient is constant across the range of
quantiles $\mathcal{T}$. The examples
above test the NFCI coefficient over $\mathcal{T} = (0, 1)$ with the discretized grid of
$\lbrace 0.1, 0.25, 0.5, 0.75, 0.9 \rbrace$. One can also specify a range and the number of grid points,
such as `arg_tau_grid = list(range = c(0, 0.25), n = 5)`.

### Monotonicity test

```r
qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, adrian, "NFCI", "ker-har",
               arg_test = "monotonicity_adjacent",
               arg_alternative = "decreasing")
```

This example tests whether the NFCI coefficient is non-decreasing across quantiles, using
adjacent differences on the default grid. To test whether the coefficient is non-increasing
across quantiles, use `arg_alternative = "increasing"`.
`"monotonicity_allpair"` tests the same null using all pairs of grid quantiles.

## Documentation

Worked examples using the Adrian, Boyarchenko and Giannone (2019) data are available through
`qrfs_example("stacked")`, `qrfs_example("sup")`, and `qrfs_example("shape")`.
See `?qrfstest_stacked`, `?qrfstest_sup`, `?qrfstest_shape` and the
vignette source in `vignettes/QRfs.Rnw`.

## Data

The data `adrian` are the Adrian, Boyarchenko and Giannone (2019) one-quarter-ahead regression
data (openICPSR 113169, CC BY 4.0).
