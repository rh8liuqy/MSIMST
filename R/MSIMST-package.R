#' The 'MSIMST' package.
#'
#' @description A Bayesian monotonic single-index mixed-effect model incorporating a multivariate skew-t likelihood with survey weights adjustments. This package includes a simulation program and the associated Gibbs sampler. The single-index function is modeled as a monotonic increasing function, with a tailored Gaussian process prior to ensure accurate estimation. Random effects are assumed to follow the canonical skew-t distribution, while residuals are modeled using the multivariate Student-t distribution. Additionally, the package provides Bayesian adjustment for survey weight information.
#'
#' @docType package
#' @name MSIMST
"_PACKAGE"
#' @useDynLib MSIMST, .registration = TRUE
#' @import Rcpp
#' @importFrom stats rbinom rnorm rpois pnorm rgamma pgamma qgamma qnorm runif model.matrix plogis
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom fields rdist
#' @importFrom parallel mclapply
#' @importFrom truncnorm rtruncnorm
#' @importFrom MASS mvrnorm
#' @importFrom Rdpack reprompt
NULL