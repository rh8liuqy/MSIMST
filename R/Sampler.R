#' The Associated Gibbs Sampler
#'
#' @param X The list of design matrix.
#' @param y The list of response values.
#' @param group_info The group information for the grouped horseshoe prior. Use 0 to represent the variables with the normal priors. Use 1,2,... to present the variables with the grouped horseshoe priors. For example, c(0,0,1,1,2,3) represents first two variables with the normal prior, third and fourth variables belong to the same group with one grouped horseshoe prior, and fifth and sixth variables belong to two different groups with two independent horseshoe prior.
#' @param beta_value The initial value for the covariates' coefficients.
#' @param beta_prior_variance The variance value of the normal prior.
#' @param beta_b_value The slope parameter.
#' @param beta_lambdasq_value The first hyperparameter associated with the grouped horseshoe prior.
#' @param beta_tausq_value The second hyperparameter associated with the grouped horseshoe prior.
#' @param xi_value The parameters associated with the single index function.
#' @param xi_lengthscale_value The first hyperparameter associated with the Gaussian process kernel.
#' @param xi_tausq_value The second hyperparameter associated with the Gaussian process kernel.
#' @param g_func_type The type of priors on the single index function. Must be one of "GP" and "BP".
#' @param dsq_value The initial value of the conditional variance of the random effects.
#' @param sigmasq_value The initial value of the conditional variance of the fixed effects.
#' @param delta_value The initial value of the skewness parameter.
#' @param nu_value The initial value of the degree of freedom. Must be larger than 2.
#' @param U_value The initial values of the latent variable U. The length of \verb{U_value} must be as the same as the number of subjects.
#' @param S_value The initial values of the latent variable S. The length of \verb{S_value} must be as the same as the number of subjects.
#' @param loglik_type The type of the log-likelihood. Must be one of "skewT","skewN", and "N".
#' @param gof_K The first hyperparameter associated with the goodness of fit test. Check \insertCite{yuan2012goodness}{MSIMST} for details.
#' @param gof_L The second hyperparameter associated with the goodness of fit test. Check \insertCite{yuan2012goodness}{MSIMST} for details.
#' @param iter_warmup The number of warm-up iterations of the Gibb samplers.
#' @param iter_sampling The number of post warm-up iterations of the Gibb samplers.
#' @param verbatim TRUE/FALSE. If verbatim is TRUE, then the updating message of the Gibbs sampler will be printed.
#' @param update An integer. For example, if \verb{update} = 10, for each 10 iteration, one udpating message of the Gibbs sampler will be printed.
#' @param incremental_output TRUE/FALSE. If \verb{incremental_output} is TRUE, an incremental output will be saved. This option should not be enabled unless users anticipate the sampling process will take longer than days.
#' @param incremental_output_filename The filename of the incremental output.
#' @param incremental_output_update An integer. For example, if \verb{incremental_output_update} = 10 then for each 10 iteration, the intermediate result will be updated once. 
#' @param n_core The number of cores will be used during the Gibbs sampler. For the Windows operating system, \verb{n_core} must be 1.
#'
#' @return A list of random quantitiles drawn from the posterior distribution using the Gibbs sampler.
#' @description
#' This is the Gibbs sampler associated with the proposed single-index mixed-effects model. This Gibbs sampler supports three different likelihoods, normal, skew-normal and skew-t likelihoods and two types of priors for the single-index funcion: the Gaussian process (GP) prior and the bernstein polynomial (BP) prior.
#' 
#' @details
#' The details of the ST-GP model can be found in the vignette. Users can access the vignette using \verb{vignette(package = "MSIMST")}.
#' 
#' @examples
#' # Set the random seed.
#' set.seed(100)
#' 
#' # Simulate the data.
#' simulated_data <- reg_simulation1(N = 50,
#'                                   ni_lambda = 8,
#'                                   beta = c(0.5,0.5,0.5),
#'                                   beta_b = 1.5,
#'                                   dsq = 0.1,
#'                                   sigmasq = 0.5,
#'                                   delta = 0.6,
#'                                   nu = 5.89)
#' 
#' y <- simulated_data$y
#' X <- simulated_data$X
#' 
#' group_info <- c(0,0,0)
#' # The number of grids (L) for approximating the single index function
#' L <- 50
#' N <- length(y)
#' GP_MCMC_output <- Gibbs_Sampler(X = X,
#'                                 y = y,
#'                                 group_info = group_info,
#'                                 beta_value = c(0.5,0.5,0.5),
#'                                 beta_prior_variance = 10,
#'                                 beta_b_value = 1.5,
#'                                 beta_lambdasq_value = 1,
#'                                 beta_tausq_value = 1,
#'                                 xi_value = abs(rnorm(n = L + 1)),
#'                                 xi_lengthscale_value = 1.0,
#'                                 xi_tausq_value = 1.0,
#'                                 g_func_type = "GP",
#'                                 dsq_value = 1,
#'                                 sigmasq_value = 1,
#'                                 delta_value = 0.6,
#'                                 nu_value = 5.89,
#'                                 U_value = abs(rnorm(N)),
#'                                 S_value = abs(rnorm(N)),
#'                                 loglik_type = "skewT",
#'                                 gof_K = 10,
#'                                 gof_L = 5,
#'                                 iter_warmup = 10,
#'                                 iter_sampling = 20,
#'                                 verbatim = TRUE,
#'                                 update = 10,
#'                                 incremental_output = FALSE,
#'                                 incremental_output_filename = NULL,
#'                                 incremental_output_update = 1e6,
#'                                 n_core = 1)
#' @export
#' 
Gibbs_Sampler <- function(X,
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
                          incremental_output_update = 1e6,
                          n_core = 1) {
  # In the Windows operating system, n_core must be 1
  if (.Platform$OS.type == "windows") {
    n_core <- 1
    print("Parallel computing is not currently support in the Windows operating system.")
    print("n_core is resetted to be 1.")
  }
  # check the incremental output filename
  if (incremental_output & is.null(incremental_output_filename)) {
    stop("The incremental_output option is enabled. \nincremental_output_filename can not be NULL !")
  }
  
  p <- ncol(X[[1]])
  N <- length(y)
  L <- length(xi_value) - 1
  u <- seq(-1,1,length.out = L + 1)
  iter_total <- iter_warmup + iter_sampling
  # check the values of loglik_type
  if (! loglik_type %in% c("skewT","skewN","N")) {
    stop('loglik_type must be one of c("skewT","skewN","N") !')
  }
  
  if (loglik_type != "skewT") {
    nu_value <- 300
  }
  if (loglik_type == "N") {
    delta_value <- 0.0
  }
  
  # ensure nu_value is larger than 2
  if (nu_value < 2) {
    print("nu_value must be larger than 2!")
    print("reset nu_value to be 2 + 1e-6.")
    nu_value <- 2 + 1e-6
  }
  
  # check value of g_func_type
  if (! g_func_type %in% c("GP","BP")) {
    stop('g_func_type must be one of c("GP","BP") !')
  }
  
  # determine the values of intercept_value and beta_a_value
  intercept_value <- 0.0
  beta_a_value <- 0.0
  
  # initialize values
  beta_value_sta <- beta_value/norm(beta_value,"2")
  b_scalar_value <- b.func(nu_value)
  if (length(U_value) != N | length(S_value) != N) {
    stop("Check the dimension of initial values of U_value and S_value !")
  }
  # U must be non-negative
  U_value <- abs(U_value)
  # S must be non-negative
  S_value <- abs(S_value)
  
  g_value <- vector(mode = "list", length = N)
  for (i in 1:N) {
    Xi <- X[[i]]
    ni <- 2*nrow(Xi)
    etai <- as.numeric(Xi %*% beta_value_sta)
    if (g_func_type == "GP") {
      phiXi <- phiX_c(Xbeta = etai, u = u, L = L)
    }
    else {
      phiXi <- Dbeta(L = L + 1, t = etai)
    }
    
    g_value[[i]] <- as.numeric(phiXi %*% xi_value)
  }
  
  # create output of the sampler
  S_output <- matrix(data = NA,
                     nrow = iter_total,
                     ncol = N,
                     byrow = TRUE)
  U_output <- matrix(data = NA,
                     nrow = iter_total,
                     ncol = N,
                     byrow = TRUE)
  B_output <- matrix(data = NA,
                     nrow = iter_total,
                     ncol = N,
                     byrow = TRUE)
  sigmasq_output <- numeric(iter_total)
  dsq_output <- numeric(iter_total)
  delta_output <- numeric(iter_total)
  intercept_output <- numeric(iter_total)
  xi_output <- matrix(data = NA,
                      nrow = iter_total,
                      ncol = L + 1,
                      byrow = TRUE)
  xi_lengthscale_output <- numeric(iter_total)
  xi_tausq_output <- numeric(iter_total)
  beta_sta_output <- matrix(data = NA,
                            nrow = iter_total,
                            ncol = p,
                            byrow = TRUE)
  beta_a_output <- numeric(iter_total)
  beta_b_output <- numeric(iter_total)
  num_group <- max(group_info)
  
  # if the horseshoe prior is used, check the value of global shrinkage
  if (num_group > 1) {
    if (beta_tausq_value < 0 | beta_tausq_value > 1) {
      warning("beta_tausq_value must be between 0 and 1.\nChange beta_tausq_value to be 0.5.")
      beta_tausq_value <- 0.5
    } 
  }
  beta_lambdasq_output <- matrix(data = NA,
                                 nrow = iter_total,
                                 ncol = max(group_info))
  beta_tausq_output <- rep(NA,iter_total)
  nu_output <- numeric(iter_total)
  gof_output <- numeric(iter_total)
  log_lik_output <- matrix(data = NA,
                           nrow = iter_total,
                           ncol = N,
                           byrow = TRUE)
  
  for (i in 1:iter_total) {
    if (verbatim & i%%update == 0) {
      print(paste0("i:",i))
    }
    
    # update S ----------------------------------------------------------------
    
    S_value <- update_S(y = y,
                        intercept_value = intercept_value,
                        beta_a_value = beta_a_value,
                        beta_b_value = beta_b_value,
                        g_value = g_value,
                        b_scalar_value = b_scalar_value,
                        delta_value = delta_value,
                        dsq_value = dsq_value,
                        sigmasq_value = sigmasq_value,
                        U_value = U_value,
                        n_core = n_core)
    
    # update U ----------------------------------------------------------------
    
    U_value <- update_U(y = y,
                        intercept_value = intercept_value,
                        beta_a_value = beta_a_value,
                        beta_b_value = beta_b_value,
                        g_value = g_value,
                        b_scalar_value = b_scalar_value,
                        delta_value = delta_value,
                        dsq_value = dsq_value,
                        sigmasq_value = sigmasq_value,
                        nu_value = nu_value,
                        S_value = S_value,
                        n_core = n_core)
    
    
    # update B ----------------------------------------------------------------
    
    B_value <- update_B(y = y,
                        intercept_value = intercept_value,
                        beta_a_value = beta_a_value,
                        beta_b_value = beta_b_value,
                        g_value = g_value,
                        b_scalar_value = b_scalar_value,
                        dsq_value = dsq_value,
                        sigmasq_value = sigmasq_value,
                        delta_value = delta_value,
                        S_value = S_value,
                        U_value = U_value)
    
    
    # update sigmasq ----------------------------------------------------------
    
    sigmasq_value <- update_sigmasq(y = y,
                                    intercept_value = intercept_value,
                                    beta_a_value = beta_a_value,
                                    beta_b_value = beta_b_value,
                                    g_value = g_value,
                                    U_value = U_value,
                                    B_value = B_value)
    
    # update dsq --------------------------------------------------------------
    
    dsq_value <- update_dsq(S_value = S_value,
                            U_value = U_value,
                            B_value = B_value,
                            b_scalar_value = b_scalar_value,
                            delta_value = delta_value)
    
    
    # update delta ------------------------------------------------------------
    
    if (loglik_type != "N") {
      delta_value <- update_delta(y = y,
                                  intercept_value = intercept_value,
                                  beta_a_value = beta_a_value,
                                  beta_b_value = beta_b_value,
                                  g_value = g_value,
                                  b_scalar_value = b_scalar_value,
                                  delta_value = delta_value,
                                  dsq_value = dsq_value,
                                  sigmasq_value = sigmasq_value,
                                  S_value = S_value,
                                  U_value = U_value,
                                  n_core = n_core)
    }
    
    # update beta_b -----------------------------------------------------------
    
    beta_b_value <- update_beta_b(y = y,
                                  intercept_value = intercept_value,
                                  beta_a_value = beta_a_value,
                                  beta_b_value = beta_b_value,
                                  g_value = g_value,
                                  delta_value = delta_value,
                                  dsq_value = dsq_value,
                                  sigmasq_value = sigmasq_value,
                                  nu_value = nu_value,
                                  S_value = S_value,
                                  U_value = U_value,
                                  B_value = B_value)
    
    # update xi ---------------------------------------------------------------
    
    if (g_func_type == "GP") {
      xi_value <- update_xi(y = y,
                            X = X,
                            intercept_value = intercept_value,
                            beta_value = beta_value,
                            beta_a_value = beta_a_value,
                            beta_b_value = beta_b_value,
                            xi_value = xi_value,
                            xi_lengthscale_value = xi_lengthscale_value,
                            xi_tausq_value = xi_tausq_value,
                            b_scalar_value = b_scalar_value,
                            delta_value = delta_value,
                            dsq_value = dsq_value,
                            sigmasq_value = sigmasq_value,
                            S_value = S_value,
                            U_value = U_value,
                            u = u,
                            L = L,
                            n_core = n_core)
      
      # update xi_lengthscale
      xi_lengthscale_value <- update_xi_lengthscale(xi_value = xi_value,
                                                    xi_lengthscale_value = xi_lengthscale_value,
                                                    xi_tausq_value = xi_tausq_value)
      
      # update xi-tausq
      xi_tausq_value <- update_xi_tausq(xi_value = xi_value,
                                        xi_lengthscale_value = xi_lengthscale_value,
                                        xi_tausq_value = xi_tausq_value)
    } else {
      xi_value <- update_xi_BP(y = y,
                               X = X,
                               intercept_value = intercept_value,
                               beta_value = beta_value,
                               beta_a_value = beta_a_value,
                               beta_b_value = beta_b_value,
                               xi_value = xi_value,
                               xi_lengthscale_value = xi_lengthscale_value,
                               xi_tausq_value = xi_tausq_value,
                               b_scalar_value = b_scalar_value,
                               delta_value = delta_value,
                               dsq_value = dsq_value,
                               sigmasq_value = sigmasq_value,
                               S_value = S_value,
                               U_value = U_value,
                               u = u,
                               n_core = n_core)
      
      # update xi_lengthscale
      xi_lengthscale_value <- update_xi_lengthscale(xi_value = xi_value,
                                                    xi_lengthscale_value = xi_lengthscale_value,
                                                    xi_tausq_value = xi_tausq_value)
      
      # update xi-tausq
      xi_tausq_value <- update_xi_tausq(xi_value = xi_value,
                                        xi_lengthscale_value = xi_lengthscale_value,
                                        xi_tausq_value = xi_tausq_value)
    }
    
    # update beta -------------------------------------------------------------
    
    if (g_func_type == "GP") {
      beta_value <- update_beta(y = y,
                                X = X,
                                group_info = group_info,
                                intercept_value = intercept_value,
                                beta_value = beta_value,
                                beta_a_value = beta_a_value,
                                beta_b_value = beta_b_value,
                                beta_lambdasq_value = beta_lambdasq_value,
                                beta_tausq_value = beta_tausq_value,
                                xi_value = xi_value,
                                dsq_value = dsq_value,
                                sigmasq_value = sigmasq_value,
                                delta_value = delta_value,
                                nu_value = nu_value,
                                u = u,
                                L = L,
                                beta_prior_variance = beta_prior_variance,
                                n_core = n_core)
      beta_value_sta <- beta_value/norm(beta_value,"2")
    } else {
      beta_value <- update_beta_BP(y = y,
                                   X = X,
                                   group_info = group_info,
                                   intercept_value = intercept_value,
                                   beta_value = beta_value,
                                   beta_a_value = beta_a_value,
                                   beta_b_value = beta_b_value,
                                   beta_lambdasq_value = beta_lambdasq_value,
                                   beta_tausq_value = beta_tausq_value,
                                   xi_value = xi_value,
                                   dsq_value = dsq_value,
                                   sigmasq_value = sigmasq_value,
                                   delta_value = delta_value,
                                   nu_value = nu_value,
                                   u = u,
                                   beta_prior_variance = beta_prior_variance,
                                   n_core = n_core)
      beta_value_sta <- beta_value/norm(beta_value,"2")
    }
    
    # update g_value - after update beta
    g_value <- mclapply(X = 1:N,
                        FUN = function(j) {
                          Xi <- X[[j]]
                          etai <- as.numeric(Xi %*% beta_value_sta)
                          if (g_func_type == "GP") {
                            phiXi <- phiX_c(Xbeta = etai, u = u, L = L)
                          } else {
                            phiXi <- Dbeta(L = L + 1, t = etai)
                          }
                          output <- as.numeric(phiXi %*% xi_value)
                          return(output)
                        },
                        mc.cores = n_core)
    
    # update horseshoe prior hyper parameters
    if (num_group > 0) {
      beta_lambdasq_value <- update_beta_lambdasq(group_info = group_info,
                                                  beta_value = beta_value,
                                                  beta_lambdasq_value = beta_lambdasq_value,
                                                  beta_tausq_value = beta_tausq_value)
      
      beta_tausq_value <- update_beta_tausq(group_info = group_info,
                                            beta_value = beta_value,
                                            beta_lambdasq_value = beta_lambdasq_value,
                                            beta_tausq_value = beta_tausq_value)
    }
    
    # update nu ---------------------------------------------------------------
    
    if (loglik_type == "skewT") {
      nu_value <- update_nu(y = y,
                            X = X,
                            intercept_value = intercept_value,
                            beta_a_value = beta_a_value,
                            beta_b_value = beta_b_value,
                            g_value = g_value,
                            dsq_value = dsq_value,
                            sigmasq_value = sigmasq_value,
                            delta_value = delta_value,
                            nu_value = nu_value,
                            n_core = n_core)
      # update b_scalar_value
      b_scalar_value <- b.func(nu_value)
    }
    else {
      b_scalar_value <- -sqrt(2/pi) # upper limit of b.func
    }
    
    
    # update goodness of fit statistics ---------------------------------------
    update_gof_value <- update_gof(y = y,
                                   g_value = g_value,
                                   intercept_value = intercept_value,
                                   beta_a_value = beta_a_value,
                                   beta_b_value = beta_b_value,
                                   B_value = B_value,
                                   U_value = U_value,
                                   delta_value = delta_value,
                                   dsq_value = dsq_value,
                                   sigmasq_value = sigmasq_value,
                                   nu_value = nu_value,
                                   b_scalar_value = b_scalar_value,
                                   gof_K = gof_K,
                                   gof_L = gof_L,
                                   n_core = n_core)
    gof_value <- update_gof_value$d
    log_lik_value <- update_gof_value$log_lik
    
    # update output values
    S_output[i,] <- S_value
    U_output[i,] <- U_value
    B_output[i,] <- B_value
    sigmasq_output[i] <- sigmasq_value
    dsq_output[i] <- dsq_value
    delta_output[i] <- delta_value
    intercept_output[i] <- intercept_value
    beta_a_output[i] <- beta_a_value
    beta_b_output[i] <- beta_b_value
    xi_output[i,] <- xi_value
    xi_lengthscale_output[i] <- xi_lengthscale_value
    xi_tausq_output[i] <- xi_tausq_value
    beta_sta_output[i,] <- beta_value_sta
    if (num_group > 0) {
      beta_lambdasq_output[i,] <- beta_lambdasq_value
      beta_tausq_output[i] <- beta_tausq_value
    }
    nu_output[i] <- nu_value
    gof_output[i] <- gof_value
    log_lik_output[i,] <- log_lik_value
    
    # save incremental output -------------------------------------------------
    
    if (incremental_output & i > 0 & i %% incremental_output_update == 0) {
      incremental_output_temp <- list(S_output = S_output[1:i,],
                                      U_output = U_output[1:i,],
                                      B_output = B_output[1:i,],
                                      sigmasq_output = sigmasq_output[1:i],
                                      dsq_output = dsq_output[1:i],
                                      delta_output = delta_output[1:i],
                                      intercept_output = intercept_output[1:i],
                                      beta_a_output = beta_a_output[1:i],
                                      beta_b_output = beta_b_output[1:i],
                                      xi_output = xi_output[1:i,],
                                      xi_lengthscale_output = xi_lengthscale_output[1:i],
                                      xi_tausq_output = xi_tausq_output[1:i],
                                      beta_output = beta_sta_output[1:i,],
                                      beta_lambdasq_output = beta_lambdasq_output[1:i,],
                                      beta_tausq_output = beta_tausq_output[1:i],
                                      nu_output = nu_output[1:i],
                                      gof_output = gof_output[1:i],
                                      log_lik_output = log_lik_output[1:i,],
                                      g_func_type = g_func_type,
                                      loglik_type = loglik_type,
                                      gof_K = gof_K,
                                      gof_L = gof_L)
      saveRDS(incremental_output_temp,incremental_output_filename)
      rm(incremental_output_temp)
    }
    
  }
  if (num_group > 0) {
    output <- list(S_output = S_output[(iter_warmup+1):iter_total,],
                   U_output = U_output[(iter_warmup+1):iter_total,],
                   B_output = B_output[(iter_warmup+1):iter_total,],
                   sigmasq_output = sigmasq_output[(iter_warmup+1):iter_total],
                   dsq_output = dsq_output[(iter_warmup+1):iter_total],
                   delta_output = delta_output[(iter_warmup+1):iter_total],
                   intercept_output = intercept_output[(iter_warmup+1):iter_total],
                   beta_a_output = beta_a_output[(iter_warmup+1):iter_total],
                   beta_b_output = beta_b_output[(iter_warmup+1):iter_total],
                   xi_output = xi_output[(iter_warmup+1):iter_total,],
                   xi_lengthscale_output = xi_lengthscale_output[(iter_warmup+1):iter_total],
                   xi_tausq_output = xi_tausq_output[(iter_warmup+1):iter_total],
                   beta_output = beta_sta_output[(iter_warmup+1):iter_total,],
                   beta_lambdasq_output = beta_lambdasq_output[(iter_warmup+1):iter_total,],
                   beta_tausq_output = beta_tausq_output[(iter_warmup+1):iter_total],
                   nu_output = nu_output[(iter_warmup+1):iter_total],
                   gof_output = gof_output[(iter_warmup+1):iter_total],
                   log_lik_output = log_lik_output[(iter_warmup+1):iter_total,])
  }
  else {
    output <- list(S_output = S_output[(iter_warmup+1):iter_total,],
                   U_output = U_output[(iter_warmup+1):iter_total,],
                   B_output = B_output[(iter_warmup+1):iter_total,],
                   sigmasq_output = sigmasq_output[(iter_warmup+1):iter_total],
                   dsq_output = dsq_output[(iter_warmup+1):iter_total],
                   delta_output = delta_output[(iter_warmup+1):iter_total],
                   intercept_output = intercept_output[(iter_warmup+1):iter_total],
                   beta_a_output = beta_a_output[(iter_warmup+1):iter_total],
                   beta_b_output = beta_b_output[(iter_warmup+1):iter_total],
                   xi_output = xi_output[(iter_warmup+1):iter_total,],
                   xi_lengthscale_output = xi_lengthscale_output[(iter_warmup+1):iter_total],
                   xi_tausq_output = xi_tausq_output[(iter_warmup+1):iter_total],
                   beta_output = beta_sta_output[(iter_warmup+1):iter_total,],
                   nu_output = nu_output[(iter_warmup+1):iter_total],
                   gof_output = gof_output[(iter_warmup+1):iter_total],
                   log_lik_output = log_lik_output[(iter_warmup+1):iter_total,])
  }
  return(output)
}
