#' The 'MSIMST' package.
#'
#' @description Incorporates a Bayesian monotonic single-index mixed-effect model with a multivariate skew-t likelihood, specifically designed to handle survey weights adjustments. Features include a simulation program and an associated Gibbs sampler for model estimation. The single-index function is constrained to be monotonic increasing, utilizing a customized Gaussian process prior for precise estimation. The model assumes random effects follow a canonical skew-t distribution, while residuals are represented by a multivariate Student-t distribution. Offers robust Bayesian adjustments to integrate survey weight information effectively.
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