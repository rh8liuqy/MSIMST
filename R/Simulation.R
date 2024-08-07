#' The Generation of the Multivariate Skew-Normal Distribution
#'
#' @param mu_vec Location parameter, n by 1.
#' @param Sigma_mat Covariance matrix, n by n.
#' @param delta_vec Skewness parameter: n by 1.
#'
#' @return  One random quantity from the multivariate skew-normal distribution.
#' @noRd
rmsn <- function(mu_vec,
                 Sigma_mat,
                 delta_vec) {
  x0 <- rnorm(n = 1, mean = 0, sd = 1)
  x1 <- rmvnorm(n = 1, mean = mu_vec, sigma = Sigma_mat, method = "svd")
  x1 <- as.numeric(x1)
  output <- delta_vec * abs(x0) + x1
  return(output)
}

#' The Generation of the Multivariate Skew-t Distribution
#'
#' @param mu_vec Location parameter, n by 1.
#' @param Sigma_mat Covariance matrix, n by n.
#' @param delta_vec Skewness parameter: n by 1.
#' @param nu The degree of freedom.
#'
#' @return One random quantity from the multivariate skew-t distribution.
#' @noRd
rmst <- function(mu_vec,
                 Sigma_mat,
                 delta_vec,
                 nu) {
  n <- length(mu_vec)
  x <- rmsn(mu_vec = rep(0,n),
            Sigma_mat = Sigma_mat,
            delta_vec = delta_vec)
  u <- rgamma(n = 1, shape = nu/2, rate = nu/2)
  output <- mu_vec + u^(-0.5)*x
  return(output)
}

#' The b Function
#'
#' @param nu The degree of freedom.
#' @description
#' The b function is defined under Equation (9) from https://arxiv.org/pdf/2109.12152 .
#' 
#'
#' @return A numeric value.
#' @noRd
b.func <- function(nu) {
  output <- -sqrt(nu/pi)*gamma(0.5*(nu-1))/gamma(0.5*nu)
  return(output)
}

#' The true single-index function used in the simulation.
#'
#' @param x The single index value.
#'
#' @return The single-index function evaluated at x.
#' @noRd
true.g <- function(x){
  y <- (x+1)/2
  output <- 5*(pnorm(y, mean=0.5, sd=0.1) - pnorm(0, mean = 0.5, sd = 0.1))
  return(output)
}

#' The Function for the Simulation Study without the Variable Selection
#'
#' @param N The number of subjects.
#' @param ni_lambda The mean of Poisson distribution.
#' @param beta A 3 by 1 vector.
#' @param beta_b The slope of PD response.
#' @param dsq A part of covariance parameter.
#' @param sigmasq A part of covariance parameter.
#' @param delta The skewness parameter.
#' @param nu The degree of freedom.
#' @return A simulated dataset with the response variable \verb{y} and the design matrix \verb{X}.
#' @description
#' This is a simply simulation study that is designed to demonstrate the correctness of the proposed Gibbs sampler, \verb{Gibbs_Sampler()}.
#' @details
#' More details of the design of this simulation study can be found in the vignette. Users can access the vignette by the command \verb{vignette(package = "MSIMST")}.
#' 
#' @examples
#' set.seed(100)
#' simulated_data <- reg_simulation1(N = 50,
#'                                   ni_lambda = 8,
#'                                   beta = c(0.5,0.5,0.5),
#'                                   beta_b = 1.5,
#'                                   dsq = 0.1,
#'                                   sigmasq = 0.5,
#'                                   delta = 0.6,
#'                                   nu = 5.89)
#' y <- simulated_data$y
#' X <- simulated_data$X
#' print(head(y))
#' print(head(X))
#' 
#' @export
reg_simulation1 <- function(N,
                            ni_lambda,
                            beta,
                            beta_b,
                            dsq,
                            sigmasq,
                            delta,
                            nu) {
  output <- list(y = vector(mode = "list", length = N),
                 X = vector(mode = "list", length = N))
  p <- length(beta)
  
  if (p != 3) {
    stop("The dimension of beta must be 3")
  }
  
  b_scalar_value <- b.func(nu)
  # ensure the norm of beta is 1
  beta <- beta / norm(beta, "2")
  for (i in 1:N) {
    # ni: the number of measurements in one subject
    ni <- rpois(n = 1, lambda = ni_lambda) + 2 
    # Xi: design matrix
    X1 <- rnorm(n = ni)
    Zi <- rbinom(n =1, size = 1, prob = 0.5)
    if (Zi == 1) {
      X2 <- rnorm(n = ni, mean = 1)
    } else {
      X2 <- rnorm(n = ni, mean = -1)
    }
    Xi <- cbind(X1,X2,rep(Zi,ni))
    output$X[[i]] <- Xi
  }
  X <- do.call(rbind,output$X)
  X.scaled <- scale(X, center = TRUE, scale = TRUE)
  weights <- max(sqrt(rowSums(X.scaled*X.scaled)))
  X.std <- (X.scaled/weights)
  row_n <- 1
  intercept <- 0.0
  beta_a <- 0.0
  for (i in 1:N) {
    Xi <- output$X[[i]]
    ni <- nrow(Xi)
    Xi_std <- X.std[row_n:(row_n + ni - 1),]
    row_n <- row_n + ni
    etai <- as.numeric(Xi_std %*% beta)
    gi <- true.g(etai)
    loci_CAL <- intercept + gi + b_scalar_value*delta
    loci_PD <- beta_a + beta_b * (intercept + gi) + b_scalar_value*delta
    loci <- c(loci_CAL,loci_PD)
    Psii_mat <- matrix(dsq,nrow = 2*ni, ncol = 2*ni) + diag(rep(sigmasq,2*ni))
    deltai_vec <- rep(delta,2*ni)
    yi <- rmst(mu_vec = loci,
               Sigma_mat = Psii_mat,
               delta_vec = deltai_vec,
               nu = nu)
    output$y[[i]] <- yi
    output$X[[i]] <- Xi_std
  }
  
  return(output)
}

#' The Function for the Simulation Study with the Variable Selection
#'
#' @param N The number of subjects.
#' @param ni_lambda The mean of Poisson distribution.
#' @param beta The covariates' coefficients. A 10 by 1 vector.
#' @param beta_b The slope of PD response.
#' @param dsq A part of covariance parameter.
#' @param sigmasq A part of covariance parameter.
#' @param delta The skewness parameter.
#' @param nu The degree of freedom.
#'
#' @return A simulated dataset with the response variable \verb{y} and the design matrix \verb{X}.
#' @description
#' This simulation study is designed to demonstrate that using the grouped horseshoe prior can successfully separate signals from noise.
#' @details
#' More details of the design of this simulation study can be found in the vignette. Users can access the vignette by the command \verb{vignette(package = "MSIMST")}.
#' @examples
#' set.seed(200)
#' simulated_data <- reg_simulation2(N = 50,
#'                                   ni_lambda = 8,
#'                                   beta = c(rep(1,6),rep(0,4)),
#'                                   beta_b = 1.5,
#'                                   dsq = 0.1,
#'                                   sigmasq = 0.5,
#'                                   delta = 0.6,
#'                                   nu = 5.89)
#' 
#' y <- simulated_data$y
#' X <- simulated_data$X
#' 
#' @export
reg_simulation2 <- function(N,
                            ni_lambda,
                            beta,
                            beta_b,
                            dsq,
                            sigmasq,
                            delta,
                            nu) {
  
  # safe guard --------------------------------------------------------------
  
  if (length(beta) != 10) {
    stop("The length of beta must be 10.")
  }  
  
  # generate samples --------------------------------------------------------
  
  output_list <- vector("list", N)
  output <- list(y = vector(mode = "list", length = N),
                 X = vector(mode = "list", length = N))
  
  for (i in 1:N) {
    ni <- rpois(n = 1, lambda = 8) + 2
    
    # Categorical Variable 1 (normal prior) nonzero
    Ca1 <- sample(x = c("A","B"), size = 1, prob = c(0.5,0.5))
    
    # Categorical Variable 2 (normal prior) nonzero
    Ca2 <- sample(x = c("A","B"), size = 1, 
                  prob = c(0.13,1-0.13)) # replicate the prevalence of diabetes in the real data
    
    # Categorical Variable 3 (horseshoe prior) nonzero
    Ca3 <- sample(x = c("A","B","C"), size = 1, 
                  prob = c(1/3,1/3,1/3))
    
    # Categorical Variable 4 (horseshoe prior) nonzero
    Ca4 <- sample(x = c("A","B"),
                  size = 1,
                  prob = c(0.5,0.5))
    
    # Continuous Variable 1 (horseshoe prior) nonzero
    if (Ca4 == "A") {
      Co1 <- rnorm(ni, mean = 1) 
    } else {
      Co1 <- rnorm(ni, mean = -1) 
    }
    
    # Categorical Variable 5 (horseshoe prior) zero
    Ca5 <- sample(x = c("A","B","C"), size = 1, 
                  prob = c(1/3,1/3,1/3))
    
    # Categorical Variable 6 (horseshoe prior) zero
    Ca6 <- sample(x = c("A","B"),
                  size = 1,
                  prob = c(0.5,0.5))
    
    # Continuous Variable 2 (horseshoe prior) zero
    if (Ca6 == "A") {
      Co2 <- rnorm(ni, mean = 1) 
    } else {
      Co2 <- rnorm(ni, mean = -1) 
    }
    
    df_temp <- data.frame(Ca1 = rep(Ca1,ni),
                          Ca2 = rep(Ca2,ni),
                          Ca3 = rep(Ca3,ni),
                          Ca4 = rep(Ca4,ni),
                          Co1 = Co1,
                          Ca5 = rep(Ca5,ni),
                          Ca6 = rep(Ca6,ni),
                          Co2 = Co2)
    
    output_list[[i]] <- df_temp
  }
  b_scalar_value <- b.func(nu)
  beta <- beta / norm(beta, "2")
  
  df_X <- do.call("rbind", output_list)
  X <- model.matrix(~ . ,df_X)
  X <- X[,-1] # drop the intercept
  
  X.scaled <- scale(X, center = TRUE, scale = TRUE)
  weights <- max(sqrt(rowSums(X.scaled*X.scaled)))
  X.std <- (X.scaled/weights)
  row_n <- 1
  
  for (i in 1:N) {
    
    ni <- nrow(output_list[[i]])
    Xi_std <- X.std[row_n:(row_n + ni - 1),]
    row_n <- row_n + ni
    etai <- as.numeric(Xi_std %*% beta)
    gi <- true.g(etai)
    loci_CAL <- gi + b_scalar_value*delta
    loci_PD <- beta_b * gi + b_scalar_value*delta
    
    # location of Y
    loci <- c(loci_CAL,loci_PD)
    
    # covariance matrix of Y
    Psii_mat <- matrix(dsq,nrow = 2*ni, ncol = 2*ni) + diag(rep(sigmasq,2*ni))
    
    # skewness of Y
    deltai_vec <- rep(delta,2*ni)
    
    y <- rmst(mu_vec = loci,
              Sigma_mat = Psii_mat,
              delta_vec = deltai_vec,
              nu = nu)
    
    output$y[[i]] <- as.numeric(y)
    output$X[[i]] <- Xi_std
  }
  
  return(output)
  
}

#' The Function for the Simulation Study with the Variable Selection and Survey Weights
#'
#' @param N The number of subjects.
#' @param ni_lambda The mean of Poisson distribution.
#' @param beta The covariates' coefficients. A 10 by 1 vector.
#' @param beta_b The slope of PD response.
#' @param dsq A part of covariance parameter.
#' @param sigmasq A part of covariance parameter.
#' @param delta The skewness parameter.
#' @param nu The degree of freedom.
#' @param muz The location parameter of the latent/selection variable.
#' @param rho The correlation parameter of the latent/selection variable.
#' @param sigmasq_z The variance parameter of the latent/selection variable.
#' @param zeta0 The intercept term inside the logistic function.
#' @param zeta1 The slope term inside the logistic function.
#'
#' @return A simulated dataset with the response variable \verb{y}, the design matrix \verb{X} and the survey weight \verb{survey_weight}.
#' @description
#' This simulation study is designed to show the effectiveness of the grouped horseshoe prior for the variable selection and the \verb{WFPBB()} function for adjusting survey weights.
#' 
#' @details
#' More details of the design of this simulation study can be found in the vignette. Users can access the vignette by the command \verb{vignette(package = "MSIMST")}.
#' @examples
#' set.seed(100)
#' output_data <- reg_simulation3(N = 1000,
#'                                ni_lambda= 8,
#'                                beta = c(rep(1,6),rep(0,4)),
#'                                beta_b = 1.5,
#'                                dsq = 0.1,
#'                                sigmasq = 0.5,
#'                                delta = 0.6,
#'                                nu = 5.89,
#'                                muz = 0,
#'                                rho = 36.0,
#'                                sigmasq_z = 0.6,
#'                                zeta0 = -1.8,
#'                                zeta1 = 0.1)
#' y <- output_data$y
#' X <- output_data$X
#' survey_weight <- output_data$survey_weight
#' 
#' @export
reg_simulation3 <- function(N,
                            ni_lambda,
                            beta,
                            beta_b,
                            dsq,
                            sigmasq,
                            delta,
                            nu,
                            muz,
                            rho,
                            sigmasq_z,
                            zeta0,
                            zeta1) {
  
  # safe guard --------------------------------------------------------------
  
  if (length(beta) != 10) {
    stop("The length of beta must be 10.")
  }  
  
  # generate samples --------------------------------------------------------
  
  while (TRUE) {
    output_list <- vector("list", N)
    output <- list(y = vector(mode = "list", length = N),
                   X = vector(mode = "list", length = N),
                   pis = rep(NA,N))
    
    for (i in 1:N) {
      ni <- rpois(n = 1, lambda = 8) + 2
      
      # Categorical Variable 1 (normal prior) nonzero
      Ca1 <- sample(x = c("A","B"), size = 1, prob = c(0.5,0.5))
      
      # Categorical Variable 2 (normal prior) nonzero
      Ca2 <- sample(x = c("A","B"), size = 1, 
                    prob = c(0.13,1-0.13)) # replicate the prevalence of diabetes in the real data
      
      # Categorical Variable 3 (horseshoe prior) nonzero
      Ca3 <- sample(x = c("A","B","C"), size = 1, 
                    prob = c(1/3,1/3,1/3))
      
      # Categorical Variable 4 (horseshoe prior) nonzero
      Ca4 <- sample(x = c("A","B"),
                    size = 1,
                    prob = c(0.5,0.5))
      
      # Continuous Variable 1 (horseshoe prior) nonzero
      if (Ca4 == "A") {
        Co1 <- rnorm(ni, mean = 1) 
      } else {
        Co1 <- rnorm(ni, mean = -1) 
      }
      
      # Categorical Variable 5 (horseshoe prior) zero
      Ca5 <- sample(x = c("A","B","C"), size = 1, 
                    prob = c(1/3,1/3,1/3))
      
      # Categorical Variable 6 (horseshoe prior) zero
      Ca6 <- sample(x = c("A","B"),
                    size = 1,
                    prob = c(0.5,0.5))
      
      # Continuous Variable 2 (horseshoe prior) zero
      if (Ca6 == "A") {
        Co2 <- rnorm(ni, mean = 1) 
      } else {
        Co2 <- rnorm(ni, mean = -1) 
      }
      
      df_temp <- data.frame(Ca1 = rep(Ca1,ni),
                            Ca2 = rep(Ca2,ni),
                            Ca3 = rep(Ca3,ni),
                            Ca4 = rep(Ca4,ni),
                            Co1 = Co1,
                            Ca5 = rep(Ca5,ni),
                            Ca6 = rep(Ca6,ni),
                            Co2 = Co2)
      
      output_list[[i]] <- df_temp
    }
    b_scalar_value <- b.func(nu)
    beta <- beta / norm(beta, "2")
    
    df_X <- do.call("rbind", output_list)
    X <- model.matrix(~ . ,df_X)
    X <- X[,-1] # drop the intercept
    
    X.scaled <- scale(X, center = TRUE, scale = TRUE)
    weights <- max(sqrt(rowSums(X.scaled*X.scaled)))
    X.std <- (X.scaled/weights)
    row_n <- 1
    
    for (i in 1:N) {
      
      ni <- nrow(output_list[[i]])
      Xi_std <- X.std[row_n:(row_n + ni - 1),]
      row_n <- row_n + ni
      etai <- as.numeric(Xi_std %*% beta)
      gi <- true.g(etai)
      loci_CAL <- gi + b_scalar_value*delta
      loci_PD <- beta_b * gi + b_scalar_value*delta
      
      # location of Y
      loci <- c(loci_CAL,loci_PD)
      
      # location of Y and Z
      loc_YZ <- c(loci,muz)
      
      # covariance matrix of Y
      Psii_mat <- matrix(dsq,nrow = 2*ni, ncol = 2*ni) + diag(rep(sigmasq,2*ni))
      cov_YZ <- matrix(0, nrow = 2*ni + 1, ncol = 2*ni + 1)
      cov_YZ[1:(2*ni),1:(2*ni)] <- Psii_mat
      cov_YZ[2*ni + 1, 1:(2*ni)] <- rho
      cov_YZ[1:(2*ni), 2*ni + 1] <- rho
      cov_YZ[2*ni + 1, 2*ni + 1] <- sigmasq_z
      
      # skewness of Y
      deltai_vec <- rep(delta,2*ni)
      
      # skewness of Y and Z
      delta_YZ <- c(deltai_vec,0)
      
      yZ <- rmst(mu_vec = loc_YZ,
                 Sigma_mat = cov_YZ,
                 delta_vec = delta_YZ,
                 nu = nu)
      
      y <- as.numeric(yZ[1:(2*ni)])
      Z <- as.numeric(yZ[2*ni + 1])
      
      prob_pi <- plogis(zeta0 + zeta1*Z)
      
      selection_var <- rbinom(n = 1, size = 1, prob = prob_pi)
      
      if (selection_var == 1) {
        
        output$y[[i]] <- y
        output$X[[i]] <- Xi_std
        output$pis[i] <- prob_pi
        
      }
    }
    
    # remove NULL values
    output$y <- output$y[sapply(output$y, function(x){ ! is.null(x)})]
    output$X <- output$X[sapply(output$X, function(x){ ! is.null(x)})]
    output$pis <- output$pis[!is.na(output$pis)]
    
    # calculation of weight
    n <- length(output$y)
    pis <- output$pis
    sum_inv_pi <- sum(1/pis)
    n <- length(output$y)
    ws <- n * (1/pis)/sum_inv_pi
    survey_weight <- ws/sum(ws) * N
    output$survey_weight <- survey_weight
    
    # remove pis
    output$pis <- NULL
    
    # ensure survey weights are no smaller than 1
    if (all(output$survey_weight >= 1)) { 
      break
      }
  }
  
  return(output)
  
}