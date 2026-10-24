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
qrfstest_stacked(ggdp_h1 ~ ggdp + NFCI, adrian,
                 "NFCI[0.1] = NFCI[0.5] = NFCI[0.9]", "os-har")
```

The data `adrian` are the Adrian, Boyarchenko and Giannone (2019) one-quarter-ahead regression
data (openICPSR 113169, CC BY 4.0).
