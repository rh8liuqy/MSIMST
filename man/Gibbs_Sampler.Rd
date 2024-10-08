% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sampler.R
\name{Gibbs_Sampler}
\alias{Gibbs_Sampler}
\title{The Associated Gibbs Sampler}
\usage{
Gibbs_Sampler(
  X,
  y,
  group_info,
  beta_value,
  beta_prior_variance,
  beta_b_value,
  beta_lambdasq_value,
  beta_tausq_value,
  xi_value,
  xi_lengthscale_value,
  xi_tausq_value,
  g_func_type,
  dsq_value,
  sigmasq_value,
  delta_value,
  nu_value,
  U_value,
  S_value,
  loglik_type,
  gof_K,
  gof_L,
  iter_warmup,
  iter_sampling,
  verbatim,
  update = 10,
  incremental_output = FALSE,
  incremental_output_filename = NULL,
  incremental_output_update = 1e+06,
  n_core = 1
)
}
\arguments{
\item{X}{The list of design matrix.}

\item{y}{The list of response values.}

\item{group_info}{The group information for the grouped horseshoe prior. Use 0 to represent the variables with the normal priors. Use 1,2,... to present the variables with the grouped horseshoe priors. For example, c(0,0,1,1,2,3) represents first two variables with the normal prior, third and fourth variables belong to the same group with one grouped horseshoe prior, and fifth and sixth variables belong to two different groups with two independent horseshoe prior.}

\item{beta_value}{The initial value for the covariates' coefficients.}

\item{beta_prior_variance}{The variance value of the normal prior.}

\item{beta_b_value}{The slope parameter.}

\item{beta_lambdasq_value}{The first hyperparameter associated with the grouped horseshoe prior.}

\item{beta_tausq_value}{The second hyperparameter associated with the grouped horseshoe prior.}

\item{xi_value}{The parameters associated with the single index function.}

\item{xi_lengthscale_value}{The first hyperparameter associated with the Gaussian process kernel.}

\item{xi_tausq_value}{The second hyperparameter associated with the Gaussian process kernel.}

\item{g_func_type}{The type of priors on the single index function. Must be one of "GP" and "BP".}

\item{dsq_value}{The initial value of the conditional variance of the random effects.}

\item{sigmasq_value}{The initial value of the conditional variance of the fixed effects.}

\item{delta_value}{The initial value of the skewness parameter.}

\item{nu_value}{The initial value of the degree of freedom. Must be larger than 2.}

\item{U_value}{The initial values of the latent variable U. The length of \verb{U_value} must be as the same as the number of subjects.}

\item{S_value}{The initial values of the latent variable S. The length of \verb{S_value} must be as the same as the number of subjects.}

\item{loglik_type}{The type of the log-likelihood. Must be one of "skewT","skewN", and "N".}

\item{gof_K}{The first hyperparameter associated with the goodness of fit test. Check \insertCite{yuan2012goodness}{MSIMST} for details.}

\item{gof_L}{The second hyperparameter associated with the goodness of fit test. Check \insertCite{yuan2012goodness}{MSIMST} for details.}

\item{iter_warmup}{The number of warm-up iterations of the Gibb samplers.}

\item{iter_sampling}{The number of post warm-up iterations of the Gibb samplers.}

\item{verbatim}{TRUE/FALSE. If verbatim is TRUE, then the updating message of the Gibbs sampler will be printed.}

\item{update}{An integer. For example, if \verb{update} = 10, for each 10 iteration, one udpating message of the Gibbs sampler will be printed.}

\item{incremental_output}{TRUE/FALSE. If \verb{incremental_output} is TRUE, an incremental output will be saved. This option should not be enabled unless users anticipate the sampling process will take longer than days.}

\item{incremental_output_filename}{The filename of the incremental output.}

\item{incremental_output_update}{An integer. For example, if \verb{incremental_output_update} = 10 then for each 10 iteration, the intermediate result will be updated once.}

\item{n_core}{The number of cores will be used during the Gibbs sampler. For the Windows operating system, \verb{n_core} must be 1.}
}
\value{
A list of random quantitiles drawn from the posterior distribution using the Gibbs sampler.
}
\description{
This is the Gibbs sampler associated with the proposed single-index mixed-effects model. This Gibbs sampler supports three different likelihoods, normal, skew-normal and skew-t likelihoods and two types of priors for the single-index funcion: the Gaussian process (GP) prior and the bernstein polynomial (BP) prior.
}
\details{
The details of the ST-GP model can be found in the vignette. Users can access the vignette using \verb{vignette(package = "MSIMST")}.
}
\examples{
# Set the random seed.
set.seed(100)

# Simulate the data.
simulated_data <- reg_simulation1(N = 50,
                                  ni_lambda = 8,
                                  beta = c(0.5,0.5,0.5),
                                  beta_b = 1.5,
                                  dsq = 0.1,
                                  sigmasq = 0.5,
                                  delta = 0.6,
                                  nu = 5.89)

y <- simulated_data$y
X <- simulated_data$X

group_info <- c(0,0,0)
# The number of grids (L) for approximating the single index function
L <- 50
N <- length(y)
GP_MCMC_output <- Gibbs_Sampler(X = X,
                                y = y,
                                group_info = group_info,
                                beta_value = c(0.5,0.5,0.5),
                                beta_prior_variance = 10,
                                beta_b_value = 1.5,
                                beta_lambdasq_value = 1,
                                beta_tausq_value = 1,
                                xi_value = abs(rnorm(n = L + 1)),
                                xi_lengthscale_value = 1.0,
                                xi_tausq_value = 1.0,
                                g_func_type = "GP",
                                dsq_value = 1,
                                sigmasq_value = 1,
                                delta_value = 0.6,
                                nu_value = 5.89,
                                U_value = abs(rnorm(N)),
                                S_value = abs(rnorm(N)),
                                loglik_type = "skewT",
                                gof_K = 10,
                                gof_L = 5,
                                iter_warmup = 10,
                                iter_sampling = 20,
                                verbatim = TRUE,
                                update = 10,
                                incremental_output = FALSE,
                                incremental_output_filename = NULL,
                                incremental_output_update = 1e6,
                                n_core = 1)
}
