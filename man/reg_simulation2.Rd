% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation.R
\name{reg_simulation2}
\alias{reg_simulation2}
\title{The Function for the Simulation Study with the Variable Selection}
\usage{
reg_simulation2(N, ni_lambda, beta, beta_b, dsq, sigmasq, delta, nu)
}
\arguments{
\item{N}{The number of subjects.}

\item{ni_lambda}{The mean of Poisson distribution.}

\item{beta}{The covariates' coefficients. A 10 by 1 vector.}

\item{beta_b}{The slope of PD response.}

\item{dsq}{A part of covariance parameter.}

\item{sigmasq}{A part of covariance parameter.}

\item{delta}{The skewness parameter.}

\item{nu}{The degree of freedom.}
}
\value{
A simulated dataset with the response variable \verb{y} and the design matrix \verb{X}.
}
\description{
This simulation study is designed to demonstrate that using the grouped horseshoe prior can successfully separate signals from noise.
}
\details{
More details of the design of this simulation study can be found in the vignette. Users can access the vignette by the command \verb{vignette(package = "MSIMST")}.
}
\examples{
set.seed(200)
simulated_data <- reg_simulation2(N = 50,
                                  ni_lambda = 8,
                                  beta = c(rep(1,6),rep(0,4)),
                                  beta_b = 1.5,
                                  dsq = 0.1,
                                  sigmasq = 0.5,
                                  delta = 0.6,
                                  nu = 5.89)

y <- simulated_data$y
X <- simulated_data$X

}
