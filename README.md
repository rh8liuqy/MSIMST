
# MSIMST

<!-- badges: start -->
<!-- badges: end -->

The goal of MSIMST is to provide a Bayesian monotonic single-index mixed-effect model incorporating a multivariate skew-t likelihood with survey weights adjustments. This package includes a simulation program and the associated Gibbs sampler. The single-index function is modeled as a monotonic increasing function, with a tailored Gaussian process prior to ensure accurate estimation. Random effects are assumed to follow the canonical skew-t distribution, while residuals are modeled using the multivariate Student-t distribution. Additionally, the package provides Bayesian adjustment for survey weight information.

## Installation

You can install the development version of MSIMST like so:

``` r
devtools::install_github(repo = "https://github.com/rh8liuqy/MSIMST")
```

## Vignette

Users can access the vignette:

``` r
library(MSIMST)
vignette("MSIMST_vignette",package = "MSIMST")
```

