# QRfs

Fixed-smoothing (fixed-*b* and fixed-*K*) HAR inference for time-series quantile regression,
following Chen, Galvao, Rho, Vogelsang and Yoon (2026), *Fixed-smoothing Uniform Inference for
Quantile Regression*.

## Installation

```r
# install.packages("remotes")
remotes::install_github("seunghwaRho/QRfs")
```

## Example

```r
library(QRfs)
```

### Uniform sup test over a quantile range

```r
sup <- qrfstest_sup(ggdp_h1 ~ ggdp + NFCI, adrian,
                    "NFCI = 0", "ker-har",
                    arg_tau_grid = list(range = c(0, 0.25), n = 5))
uniform_CI_band(sup)
```

`qrfstest_sup()` performs a uniform-in-$\tau$ hypothesis test over a range of quantiles $\mathcal{T}$
with the sup statistic. The example above tests $H_0: \beta_{NFCI}(\tau) = 0$ for all
$\tau \in (0, 0.25]$, discretized at five grid points. When `arg_tau_grid` is not specified,
$\mathcal{T} = (0, 1)$ discretized at $\lbrace 0.1, 0.3, 0.5, 0.7, 0.9 \rbrace$. `uniform_CI_band()`
provides the uniform confidence band.
`arg_har_method = "ker-har"` chooses fixed-*b* inference with the Bartlett kernel, `"os-har"`
fixed-*K* inference with the DCT-II orthonormal series.

### Stack-Wald test at specified quantiles

```r
qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, adrian,
                 "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", "os-har")
```

`qrfstest_stacked()` performs a test of joint linear restrictions at specified quantiles
$\bar{\tau}_m$ with the stack-Wald statistic. The example above tests
$H_0: \beta_{NFCI}(0.1) = \beta_{NFCI}(0.5) = \beta_{NFCI}(0.9)$, two equalities at
$\bar{\tau}_m = \lbrace 0.1, 0.5, 0.9 \rbrace$.

### Homogeneity test

```r
qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, adrian, "NFCI", "ker-har",
               arg_test = "homogeneity")
qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, adrian, "NFCI", "ker-har",
               arg_test = "homogeneity",
               arg_tau_grid = c(0.1, 0.25, 0.5, 0.75, 0.9))
qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, adrian, "NFCI", "ker-har",
               arg_test = "homogeneity",
               arg_tau_grid = list(range = c(0.1, 0.9), n = 5))
```

`qrfstest_shape()` with `arg_test = "homogeneity"` performs the homogeneity test that one
coefficient is constant across a range of quantiles $\mathcal{T}$,
$H_0: \beta_j(\tau_1) = \beta_j(\tau_2)$ for all $\tau_1 \ne \tau_2 \in \mathcal{T}$. The examples
above test the NFCI coefficient over $\mathcal{T} = (0, 1)$, on the grid
$\lbrace 0.1, 0.25, 0.5, 0.75, 0.9 \rbrace$, and over $[0.1, 0.9]$ with five grid points. When
`arg_tau_grid` is not specified, $\mathcal{T} = (0, 1)$ discretized at
$\lbrace 0.1, 0.3, 0.5, 0.7, 0.9 \rbrace$.

### Monotonicity test

```r
qrfstest_shape(ggdp_h1 ~ ggdp + NFCI, adrian, "NFCI", "ker-har",
               arg_test = "monotonicity_adjacent",
               arg_alternative = "decreasing")
```

`qrfstest_shape()` with `arg_test = "monotonicity_adjacent"` performs the monotonicity test of
$H_0: \beta_{NFCI}(\tau)$ is non-decreasing in $\tau \in (0, 1)$ over adjacent quantiles.
`"monotonicity_allpair"` tests it over every pair of quantiles. `arg_alternative = "increasing"`
reverses the null inequality. The grid is chosen as for the homogeneity test.

## Documentation

`qrfs_example("stacked")`, `"sup"` and `"shape"` run the worked examples of the
paper's application. See `?qrfstest_stacked`, `?qrfstest_sup`, `?qrfstest_shape` and the
vignette source in `vignettes/QRfs.Rnw`.

## Data

The data `adrian` are the Adrian, Boyarchenko and Giannone (2019) one-quarter-ahead regression
data (openICPSR 113169, CC BY 4.0).
