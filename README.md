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

### Stack-Wald test at specified quantiles

```r
qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, adrian,
                 "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", "os-har")
```

`qrfstest_stacked()` performs the hypothesis test of
$H_0: \beta_{NFCI}(0.1) = \beta_{NFCI}(0.5) = \beta_{NFCI}(0.9)$ with the stack-Wald statistic.
`"ker-har"` selects fixed-*b* inference with the Bartlett kernel, `"os-har"` fixed-*K* inference
with the DCT-II orthonormal basis functions.

### Uniform test over a quantile range

```r
sup <- qrfstest_sup(ggdp_h1 ~ ggdp + NFCI, adrian,
                    "NFCI = 0", "ker-har",
                    arg_tau_grid = list(range = c(0, 0.25), n = 5))
uniform_CI_band(sup)
```

`qrfstest_sup()` performs uniform inference for $H_0: \beta_{NFCI}(\tau) = 0$ for all
$\tau \in \mathcal{T}$ with the sup statistic. When `arg_tau_grid` is not specified, $\mathcal{T} = (0, 1)$
discretized at $\lbrace 0.1, 0.3, 0.5, 0.7, 0.9 \rbrace$. The example above sets $\mathcal{T} = (0, 0.25]$
with five grid points. `uniform_CI_band()` provides the uniform confidence band.

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

`qrfstest_shape()` with `arg_test = "homogeneity"` performs the homogeneity test of
$H_0: \beta_{NFCI}(\tau)$ is constant in $\tau \in (0, 1)$. When `arg_tau_grid` is not specified,
the default grid is $\lbrace 0.1, 0.3, 0.5, 0.7, 0.9 \rbrace$ over $(0, 1)$. Supply a vector such as
`c(0.1, 0.25, 0.5, 0.75, 0.9)` to specify the grid directly, or `list(range = c(0.1, 0.9), n = 5)`
to generate a grid over a specified range.

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
