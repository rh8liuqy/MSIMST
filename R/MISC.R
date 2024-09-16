# Functions inside the "MISC.R" file are not exported to general users ----


# The Elliptical Slice Sampling -------------------------------------------

ESS <- function(loglikelihood,
                f.current,
                mu,
                Sigma) {
  p <- length(f.current)
  nu <- mvrnorm(n = 1,
                mu = mu,
                Sigma = Sigma)
  nu <- as.numeric(nu)
  u <- runif(n = 1,
             min = 0,
             max = 1)
  logy <- loglikelihood(f.current) + log(u)
  theta <- runif(n = 1,
                 min = 0,
                 max = 2*pi)
  theta.min <- theta - 2*pi
  theta.max <- theta
  condition <- TRUE
  while (condition) {
    f.new <- (f.current-mu)*cos(theta) + (nu-mu)*sin(theta) + mu
    if (loglikelihood(f.new) > logy) {
      condition <- FALSE
    }
    else {
      if (theta < 0) {
        theta.min <- theta
      }
      else {
        theta.max <- theta
      }
      theta <- runif(n = 1,
                     min = theta.min,
                     max = theta.max)
    }
  }
  return(f.new)
}

scaling_inv_matrix <- function(X_mat, 
                               beta_b_value) {
  n <- nrow(X_mat)
  X_mat[1:(n/2),(n/2+1):n] <- beta_b_value * X_mat[1:(n/2),(n/2+1):n]
  X_mat[(n/2+1):n,1:(n/2)] <- beta_b_value * X_mat[(n/2+1):n,1:(n/2)]
  X_mat[(n/2+1):n,(n/2+1):n] <- (beta_b_value^2) * X_mat[(n/2+1):n,(n/2+1):n]
  return(X_mat)
}

update_S <- function(y,
                     intercept_value,
                     beta_a_value,
                     beta_b_value,
                     g_value,
                     b_scalar_value,
                     delta_value,
                     sigmasq_value,
                     dsq_value,
                     U_value,
                     n_core) {
  N <- length(y)
  loc_vec <- numeric(N)
  scale_vec <- numeric(N)
  output.temp <- mclapply(X = 1:N,
                          FUN = function(i) {
                            yi <- y[[i]]
                            ni <- length(yi)
                            gi_value <- g_value[[i]]
                            Ui <- U_value[i]
                            phii_mat_inv_value <- cs_inv_c(dsq = dsq_value,
                                                           sigmasq = sigmasq_value,
                                                           n = ni)
                            loc1i <- intercept_value + gi_value
                            loc2i <- beta_a_value + beta_b_value*loc1i
                            loci <- c(loc1i,loc2i)
                            ystari <- yi - loci - b_scalar_value*delta_value
                            Sigmai_inv_value <- Ui*phii_mat_inv_value
                            Sigmai_inv_rowSum <- rowSums(Sigmai_inv_value)
                            p1 <- delta_value * sum(ystari * Sigmai_inv_rowSum)
                            p2 <- (delta_value^2) * sum(Sigmai_inv_value) + Ui
                            loci <- p1/p2
                            scalei <- sqrt(1/p2)
                            output <- list(loci = loci,
                                           scalei = scalei)
                          },
                          mc.cores = n_core)
  
  for (i in 1:N) {
    loc_vec[i] <- output.temp[[i]]$loci
    scale_vec[i] <- output.temp[[i]]$scalei
  }
  
  output <- rtruncnorm(n = N,
                       a = rep(0,N),
                       b = rep(Inf,N),
                       mean = loc_vec,
                       sd = scale_vec)
  return(output)
}

update_U <- function(y,
                     intercept_value,
                     beta_a_value,
                     beta_b_value,
                     g_value,
                     b_scalar_value,
                     delta_value,
                     dsq_value,
                     sigmasq_value,
                     nu_value,
                     S_value,
                     n_core) {
  N <- length(y)
  shape_vec <- numeric(N)
  rate_vec <- numeric(N)
  
  output.temp <- mclapply(X = 1:N,
                          FUN = function(i) {
                            yi <- y[[i]]
                            ni <- length(yi)
                            gi_value <- g_value[[i]]
                            Si <- S_value[i]
                            phii_mat_inv_value <- cs_inv_c(dsq = dsq_value,
                                                           sigmasq = sigmasq_value,
                                                           n = ni)
                            loc1i <- intercept_value + gi_value
                            loc2i <- beta_a_value + beta_b_value*loc1i
                            loci <- c(loc1i,loc2i)
                            ystari <- yi - loci - b_scalar_value*delta_value - delta_value*Si
                            shape_veci <- 0.5*(1 + ni + nu_value)
                            rate_veci <- 0.5*(as.numeric(t(ystari) %*% phii_mat_inv_value %*% ystari) + Si^2 + nu_value)
                            output <- list(shape_veci = shape_veci,
                                           rate_veci = rate_veci)
                            return(output)
                          },
                          mc.cores = n_core)
  
  for (i in 1:N) {
    shape_vec[i] <- output.temp[[i]]$shape_veci
    rate_vec[i] <- output.temp[[i]]$rate_veci
  }
  
  output <- rgamma(n = N,shape = shape_vec,rate = rate_vec)
  return(output)
}

update_B <- function(y,
                     intercept_value,
                     beta_a_value,
                     beta_b_value,
                     g_value,
                     b_scalar_value,
                     dsq_value,
                     sigmasq_value,
                     delta_value,
                     S_value,
                     U_value) {
  N <- length(y)
  output <- numeric(N)
  for (i in 1:N) {
    yi <- y[[i]]
    ni <- length(yi)
    gi_value <- g_value[[i]]
    Si <- S_value[i]
    Ui <- U_value[i]
    loc1i <- intercept_value + gi_value
    loc2i <- beta_a_value + beta_b_value*loc1i
    loci <- c(loc1i,loc2i)
    ystari <- yi - loci
    thetai <- delta_value * (b_scalar_value + Si)
    p1 <- (sum(ystari)/sigmasq_value + thetai/dsq_value)/(ni/sigmasq_value + 1/dsq_value)
    p2 <- 1/(Ui*(ni/sigmasq_value + 1/dsq_value))
    value <- rnorm(n = 1, mean = p1, sd = sqrt(p2))
    output[i] <- value
  }
  return(output)
}

update_sigmasq <- function(y,
                           intercept_value,
                           beta_a_value,
                           beta_b_value,
                           g_value,
                           U_value,
                           B_value) {
  N <- length(y)
  # inverse gamma prior with shape and scale
  shape <- 5
  scale <- 5
  p1 <- 0
  p2 <- 0
  for (i in 1:N) {
    yi <- y[[i]]
    ni <- length(yi)
    gi_value <- g_value[[i]]
    Ui <- U_value[i]
    Bi <- B_value[i]
    loc1i <- intercept_value + gi_value
    loc2i <- beta_a_value + beta_b_value*loc1i
    loci <- c(loc1i,loc2i)
    ystari <- yi - loci - Bi
    p1 <- p1 + 0.5*ni
    p2 <- p2 + 0.5*sum(ystari^2)*Ui
  }
  p1 <- p1 + shape
  p2 <- p2 + scale
  output <- 1/rgamma(n = 1, shape = p1, rate = p2)
  return(output)
}

update_dsq <- function(S_value,
                       U_value,
                       B_value,
                       b_scalar_value,
                       delta_value) {
  N <- length(S_value)
  # inverse gamma prior with shape and scale
  shape <- 5
  scale <- 5
  p1 <- 0
  for (i in 1:N) {
    Si <- S_value[i]
    Ui <- U_value[i]
    Bi <- B_value[i]
    Bstari <- Bi - delta_value*(b_scalar_value + Si)
    p1 <- p1 + 0.5 * (Bstari^2) * Ui
  }
  output <- 1/rgamma(n = 1, shape = N/2 + shape, rate = p1 + scale)
  return(output)
}

update_delta <- function(y,
                         intercept_value,
                         beta_a_value,
                         beta_b_value,
                         g_value,
                         b_scalar_value,
                         delta_value,
                         dsq_value,
                         sigmasq_value,
                         S_value,
                         U_value,
                         n_core) {
  N <- length(y)
  p1 <- 0
  p2 <- 0
  
  output.temp <- mclapply(X = 1:N,
                          FUN = function(i) {
                            yi <- y[[i]]
                            ni <- length(yi)
                            gi_value <- g_value[[i]]
                            Si <- S_value[i]
                            Ui <- U_value[i]
                            phii_mat_inv_value <- cs_inv_c(dsq = dsq_value,
                                                           sigmasq = sigmasq_value,
                                                           n = ni)
                            loc1i <- intercept_value + gi_value
                            loc2i <- beta_a_value + beta_b_value*loc1i
                            loci <- c(loc1i,loc2i)
                            ystari <- yi - loci
                            Vi <- (b_scalar_value + Si) * matrix(data = 1,nrow = ni,ncol = 1)
                            Sigmai_inv_value <- Ui*phii_mat_inv_value
                            Sigmai_inv_V_value <- Sigmai_inv_value %*% Vi
                            p1 <- as.numeric(t(ystari) %*% Sigmai_inv_V_value)
                            p2 <- as.numeric(t(Vi) %*% Sigmai_inv_V_value)
                            output <- list(p1 = p1,
                                           p2 = p2)
                            return(output)
                          },
                          mc.cores = n_core)
  
  for (i in 1:N) {
    p1 <- p1 + output.temp[[i]]$p1
    p2 <- p2 + output.temp[[i]]$p2
  }
  
  p2 <- 1/1000 + p2 #normal prior with mean 0 and variance 1000
  output <- rnorm(n = 1,
                  mean = p1/p2,
                  sd = sqrt(1/p2))
  return(output)
}

update_intercept <- function(y,
                             intercept_value,
                             beta_a_value,
                             beta_b_value,
                             g_value,
                             delta_value,
                             dsq_value,
                             sigmasq_value,
                             nu_value,
                             n_core) {
  logpdf <- function(intercept) {
    output <- mst_reg_lpdf_parallel(y_vec_list = y,
                                    intercept = intercept,
                                    beta_a = beta_a_value,
                                    beta_b = beta_b_value,
                                    g_value_list = g_value,
                                    dsq = dsq_value,
                                    sigmasq = sigmasq_value,
                                    delta = delta_value,
                                    nu = nu_value,
                                    n_core = n_core)
  }
  output <- ESS(loglikelihood = logpdf,
                f.current = intercept_value,
                mu = 0.0,
                Sigma = sqrt(1000)^2) # the prior variance is 1000
  return(output)
}

update_beta_a <- function(y,
                          intercept_value,
                          beta_a_value,
                          beta_b_value,
                          g_value,
                          delta_value,
                          dsq_value,
                          sigmasq_value,
                          nu_value,
                          n_core) {
  logpdf <- function(beta_a) {
    output <- mst_reg_lpdf_parallel(y_vec_list = y,
                                    intercept = intercept_value,
                                    beta_a = beta_a,
                                    beta_b = beta_b_value,
                                    g_value_list = g_value,
                                    dsq = dsq_value,
                                    sigmasq = sigmasq_value,
                                    delta = delta_value,
                                    nu = nu_value,
                                    n_core = n_core)
  }
  output <- ESS(loglikelihood = logpdf,
                f.current = beta_a_value,
                mu = 0.0,
                Sigma = sqrt(1000)^2) # the prior variance is 1000
  return(output)
}

update_beta_b <- function(y,
                          intercept_value,
                          beta_a_value,
                          beta_b_value,
                          g_value,
                          delta_value,
                          dsq_value,
                          sigmasq_value,
                          nu_value,
                          S_value,
                          U_value,
                          B_value) {
  N <- length(y)
  p1 <- 0
  p2 <- 0
  for (i in 1:N) {
    yi <- y[[i]]
    ni <- length(yi)
    gi_value <- g_value[[i]]
    Si <- S_value[i]
    Ui <- U_value[i]
    bi <- B_value[i]
    y2i <- yi[(ni/2+1):ni]
    y2istar <- y2i - bi - beta_a_value
    Sigmai_inv <- diag(x = Ui/sigmasq_value, nrow = ni/2, ncol = ni/2)
    Xi <- intercept_value + gi_value
    tXi <- t(Xi)
    tXi_Sigmai_inv <- tXi%*%Sigmai_inv
    p1 <- p1 + as.numeric(tXi_Sigmai_inv%*%y2istar)
    p2 <- p2 + as.numeric(tXi_Sigmai_inv%*%Xi)
  }
  # normal prior with mean 0 and variance 1000
  p2 <- p2 + 1/1000
  output <- rnorm(n = 1, mean = p1/p2, sd = sqrt(1/p2)) 
  return(output)
}

# Matern kernel with smoothness nu and length-scale l:
MK <- function(x, 
               y ,
               lengthscale, 
               smoothness){
  ifelse(abs(x-y)>0, 
         (sqrt(2*smoothness)*abs(x-y)/lengthscale)^smoothness/(2^(smoothness-1)*gamma(smoothness))*besselK(x=abs(x-y)*sqrt(2*smoothness)/lengthscale, nu=smoothness), 
         1.0)
}

# Covariance matrix
covmat <- function(knot,
                   smoothness,
                   lengthscale){
  return(MK(rdist(knot),
            0,
            lengthscale,
            smoothness))
}

update_xi <- function(y,
                      X,
                      intercept_value,
                      beta_value,
                      beta_a_value,
                      beta_b_value,
                      xi_value,
                      xi_lengthscale_value,
                      xi_tausq_value,
                      b_scalar_value,
                      delta_value,
                      dsq_value,
                      sigmasq_value,
                      S_value,
                      U_value,
                      u,
                      L,
                      n_core) {
  N <- length(y)
  beta_value_sta <- beta_value/norm(beta_value,"2")
  
  output_temp <- mclapply(1:N,
                          FUN = function(i) {
                            yi <- y[[i]]
                            Xi <- X[[i]]
                            ni <- length(yi)
                            Si <- S_value[i]
                            Ui <- U_value[i]
                            phii_mat_inv_value <- cs_inv_c(dsq = dsq_value,
                                                           sigmasq = sigmasq_value,
                                                           n = ni)
                            
                            # create design matrix, phiXi
                            etai <- as.numeric(Xi %*% beta_value_sta)
                            phiXi <- phiX_c(Xbeta = etai,
                                            u = u,
                                            L = L)
                            phiXi <- rbind(phiXi,phiXi)
                            tphiXi <- t(phiXi)
                            y1i <- yi[1:(ni/2)]
                            y2i <- yi[(ni/2+1):ni]
                            y1i <- y1i - intercept_value - b_scalar_value * delta_value - delta_value * Si
                            y2i <- (y2i - intercept_value*beta_b_value - beta_a_value
                                    - b_scalar_value * delta_value - delta_value * Si)/beta_b_value
                            ystari <- c(y1i,y2i)
                            Sigmai_inv <- Ui*phii_mat_inv_value
                            Sigmai_inv <- scaling_inv_matrix(X_mat = Sigmai_inv,
                                                             beta_b_value = beta_b_value)
                            tphiXi_Sigmai_inv <- tphiXi %*% Sigmai_inv
                            tphiXi_Sigmai_inv_phiXi <- tphiXi_Sigmai_inv %*% phiXi
                            p1 <- tphiXi_Sigmai_inv_phiXi
                            p2 <- tphiXi_Sigmai_inv %*% ystari
                            output <- list(p1 = p1,
                                           p2 = p2)
                            return(output)
                          },
                          mc.cores = n_core)
  
  p1 <- matrix(data = 0, nrow = nrow(output_temp[[1]]$p1), ncol = ncol(output_temp[[1]]$p1))
  p2 <- matrix(data = 0, nrow = nrow(output_temp[[1]]$p2), ncol = ncol(output_temp[[1]]$p2))
  
  for (i in 1:N) {
    p1 <- p1 + output_temp[[i]]$p1
    p2 <- p2 + output_temp[[i]]$p2
  }
  
  # prior covariance matrix
  my_knots <- seq(-1,1,length.out=L+1)
  Kmat <- covmat(knot=my_knots,
                 smoothness=1.5,
                 lengthscale=xi_lengthscale_value)
  Kmat <- xi_tausq_value * Kmat
  Kat_inv <- solve(Kmat) 
  p1 <- p1 + Kat_inv
  p1_inv <- solve(p1)
  mean_vec <- as.numeric(p1_inv %*% p2)
  output <- rtmvnormHMC(n = 1,
                        mu = mean_vec,
                        Sigma = p1_inv,
                        x_init = xi_value,
                        ff = diag(1,L + 1),
                        gg = rep(0,L + 1),
                        n_burn = 0)
  output <- as.numeric(output)
  return(output)
}

# the "prior" information for xi given hyperparameters
log_pdf_xi <- function(xi,
                       lengthscale, 
                       tausq){
  L <- length(xi) - 1
  my_knots <- seq(-1,1,length.out=L+1)
  Kmat <- covmat(knot=my_knots,
                 smoothness=1.5,
                 lengthscale=lengthscale)
  output <- dmvnorm(x = xi,
                    mean = rep(0, L + 1),
                    sigma =  tausq*Kmat,
                    log = TRUE)
  return(output)
}

update_xi_lengthscale <- function(xi_value,
                                  xi_lengthscale_value,
                                  xi_tausq_value) {
  logpdf <- function(x) {
    output <- log_pdf_xi(xi = xi_value,
                         lengthscale = exp(x), 
                         tausq = xi_tausq_value)
    return(output)
  }
  output <- ESS(loglikelihood = logpdf,
                f.current = log(xi_lengthscale_value),
                mu = 0.0,
                Sigma = 1.0^2)
  output <- exp(output)
  return(output)
}

update_xi_tausq <- function(xi_value,
                            xi_lengthscale_value,
                            xi_tausq_value) {
  logpdf <- function(x) {
    output <- log_pdf_xi(xi = xi_value,
                         lengthscale = xi_lengthscale_value, 
                         tausq = exp(x))
    return(output)
  }
  output <- ESS(loglikelihood = logpdf,
                f.current = log(xi_tausq_value),
                mu = 0.0,
                Sigma = 1.0^2)
  output <- exp(output)
  return(output)
}

update_xi_BP <- function(y,
                         X,
                         intercept_value,
                         beta_value,
                         beta_a_value,
                         beta_b_value,
                         xi_value,
                         xi_lengthscale_value,
                         xi_tausq_value,
                         b_scalar_value,
                         delta_value,
                         dsq_value,
                         sigmasq_value,
                         S_value,
                         U_value,
                         u,
                         n_core) {
  N <- length(y)
  beta_value_sta <- beta_value/norm(beta_value,"2")
  L <- length(xi_value) - 1
  p1 <- matrix(data = 0, nrow = L + 1, ncol = L + 1)
  p2 <- matrix(data = 0, nrow = L + 1, ncol = 1)
  
  output.temp <- mclapply(X = 1:N,
                          FUN = function(i) {
                            yi <- y[[i]]
                            Xi <- X[[i]]
                            ni <- length(yi)
                            Si <- S_value[i]
                            Ui <- U_value[i]
                            phii_mat_inv_value <- cs_inv_c(dsq = dsq_value,
                                                           sigmasq = sigmasq_value,
                                                           n = ni)
                            
                            # create design matrix, phiXi
                            etai <- as.numeric(Xi %*% beta_value_sta)
                            phiXi <- Dbeta(L = L + 1, t = etai)
                            phiXi <- rbind(phiXi,phiXi)
                            tphiXi <- t(phiXi)
                            y1i <- yi[1:(ni/2)]
                            y2i <- yi[(ni/2+1):ni]
                            y1i <- y1i - intercept_value - b_scalar_value * delta_value - delta_value * Si
                            y2i <- (y2i - intercept_value*beta_b_value - beta_a_value
                                    - b_scalar_value * delta_value - delta_value * Si)/beta_b_value
                            ystari <- c(y1i,y2i)
                            Sigmai_inv <- Ui*phii_mat_inv_value
                            Sigmai_inv <- scaling_inv_matrix(X_mat = Sigmai_inv,
                                                             beta_b_value = beta_b_value)
                            tphiXi_Sigmai_inv <- tphiXi %*% Sigmai_inv
                            tphiXi_Sigmai_inv_phiXi <- tphiXi_Sigmai_inv %*% phiXi
                            p1 <- tphiXi_Sigmai_inv_phiXi
                            p2 <- tphiXi_Sigmai_inv %*% ystari
                            output <- list(p1 = p1, p2 = p2)
                            return(output)
                          },
                          mc.cores = n_core)
  
  for (i in 1:N) {
    p1 <- p1 + output.temp[[i]]$p1
    p2 <- p2 + output.temp[[i]]$p2
  }
  
  # prior covariance matrix
  my_knots <- seq(-1,1,length.out=L+1)
  Kmat <- covmat(knot=my_knots,
                 smoothness=1.5,
                 lengthscale=xi_lengthscale_value)
  Kmat <- xi_tausq_value * Kmat
  Kat_inv <- solve(Kmat) 
  p1 <- p1 + Kat_inv
  p1_inv <- solve(p1)
  mean_vec <- as.numeric(p1_inv %*% p2)
  output <- rtmvnormHMC(n = 1,
                        mu = mean_vec,
                        Sigma = p1_inv,
                        x_init = xi_value,
                        ff = diag(1,L + 1),
                        gg = rep(0,L + 1),
                        n_burn = 0)
  output <- as.numeric(output)
  return(output)
}

mst_reg_lpdf_parallel <- function(y_vec_list,
                                  intercept,
                                  beta_a,
                                  beta_b,
                                  g_value_list,
                                  dsq,
                                  sigmasq,
                                  delta,
                                  nu,
                                  n_core) {
  
  # define the number of elements in a list
  N_list <- length(y_vec_list)
  output <- mclapply(1:N_list,
                     FUN = function(x) {
                       x_vec <- y_vec_list[[x]]
                       g_value <- g_value_list[[x]]
                       loc1 <- intercept + g_value
                       loc2 <- beta_a + beta_b*loc1
                       loc <- c(loc1,loc2)
                       output <- mst_lpdf_c(x_vec, loc, dsq, sigmasq, delta, nu)
                       return(output)
                     },
                     mc.cores = n_core)
  output <- Reduce("+",output)
  return(output)
}

update_beta <- function(y,
                        X,
                        group_info,
                        intercept_value,
                        beta_value,
                        beta_a_value,
                        beta_b_value,
                        beta_lambdasq_value,
                        beta_tausq_value,
                        xi_value,
                        dsq_value,
                        sigmasq_value,
                        delta_value,
                        nu_value,
                        u,
                        L,
                        beta_prior_variance,
                        n_core) {
  N <- length(y)
  p <- length(beta_value)
  logpdf <- function(beta) {
    beta_sta <- beta/norm(beta,"2")
    
    g_value_list <- mclapply(1:N,
                             FUN = function(i) {
                               Xi <- X[[i]]
                               etai <- Xi %*% beta_sta
                               phiXi <- phiX_c(etai,u,L)
                               output <- as.numeric(phiXi %*% xi_value)
                               return(output)
                             },
                             mc.cores = n_core)
    
    output <- mst_reg_lpdf_parallel(y_vec_list = y,
                                    intercept = intercept_value,
                                    beta_a = beta_a_value,
                                    beta_b = beta_b_value,
                                    g_value_list = g_value_list,
                                    dsq = dsq_value,
                                    sigmasq = sigmasq_value,
                                    delta = delta_value,
                                    nu = nu_value,
                                    n_core = n_core)
    return(output)
  }
  diag_elements <- numeric(p)
  # index of variables that must exit in the model
  index <- group_info == 0
  diag_elements[index] <- (sqrt(beta_prior_variance))^2
  # detect the number of groups for the variable selection
  num_group <- max(group_info)
  # if the variable selection is requested
  if (num_group > 0) {
    for (i in 1:num_group) {
      index <- group_info == i
      diag_elements[index] <- (beta_lambdasq_value[i])*(beta_tausq_value)
    }
  }
  Sigma <- diag(diag_elements)
  
  output <- ESS(loglikelihood = logpdf,
                f.current = beta_value,
                mu = rep(0,p),
                Sigma = Sigma)
  return(output)
}

update_beta_BP <- function(y,
                           X,
                           group_info,
                           intercept_value,
                           beta_value,
                           beta_a_value,
                           beta_b_value,
                           beta_lambdasq_value,
                           beta_tausq_value,
                           xi_value,
                           dsq_value,
                           sigmasq_value,
                           delta_value,
                           nu_value,
                           u,
                           beta_prior_variance,
                           n_core) {
  N <- length(y)
  p <- length(beta_value)
  L <- length(xi_value)
  g_value_list <- vector(mode = "list", length = N)
  logpdf <- function(beta) {
    beta_sta <- beta/norm(beta,"2")
    
    g_value_list <- mclapply(X = 1:N,
                             FUN = function(i) {
                               Xi <- X[[i]]
                               etai <- as.numeric(Xi %*% beta_sta)
                               Dmat <- Dbeta(L = L, t = etai)
                               g_value <- as.numeric(Dmat %*% xi_value)
                               return(g_value)
                             },
                             mc.cores = n_core)
    
    output <- mst_reg_lpdf_parallel(y_vec_list = y,
                                    intercept = intercept_value,
                                    beta_a = beta_a_value,
                                    beta_b = beta_b_value,
                                    g_value_list = g_value_list,
                                    dsq = dsq_value,
                                    sigmasq = sigmasq_value,
                                    delta = delta_value,
                                    nu = nu_value,
                                    n_core = n_core)
    return(output)
  }
  diag_elements <- numeric(p)
  # index of variables that must exit in the model
  index <- group_info == 0
  diag_elements[index] <- (sqrt(beta_prior_variance))^2
  # detect the number of groups for the variable selection
  num_group <- max(group_info)
  # if the variable selection is requested
  if (num_group > 0) {
    for (i in 1:num_group) {
      index <- group_info == i
      diag_elements[index] <- (beta_lambdasq_value[i])*(beta_tausq_value)
    }
  }
  Sigma <- diag(diag_elements)
  
  output <- ESS(loglikelihood = logpdf,
                f.current = beta_value,
                mu = rep(0,p),
                Sigma = Sigma)
  return(output)
}

# random number generator for (a,b) truncated gamma distribution
# with mean = shape/rate
# n: sample size
# shape: shape parameter
# rate: rate parameter
# a: lower limit of the distribution
# b: upper limit of the distribution
rtga <- function(n,shape,rate,a,b) {
  U_min <- pgamma(a, shape = shape, rate = rate)
  U_max <- pgamma(b, shape = shape, rate = rate)
  U <- runif(n = n,
             min = U_min,
             max = U_max)
  U[U>1-1e-8] <- 1-1e-8 # for numerical stability
  output <- qgamma(U, shape = shape, rate = rate)
  return(output)
}

update_beta_lambdasq <- function(group_info,
                                 beta_value,
                                 beta_lambdasq_value,
                                 beta_tausq_value) {
  num_group <- max(group_info)
  output <- numeric(num_group)
  for (i in 1:num_group) {
    index <- group_info == i
    betai <- beta_value[index]
    beta_lambdasqi <- beta_lambdasq_value[i]
    etai <- 1/beta_lambdasqi
    p <- length(betai)
    U <- runif(n = 1,
               min = 0,
               max = 1/(1+etai))
    etai <- rtga(n = 1,
                 shape = 0.5*p+0.5,
                 rate = sum(betai^2)/(2*beta_tausq_value),
                 a = 0.0,
                 b = (1-U)/U)
    lambdasqi <- 1/etai
    output[i] <- lambdasqi
  }
  return(output)
}

update_beta_tausq <- function(group_info,
                              beta_value,
                              beta_lambdasq_value,
                              beta_tausq_value) {
  # define the betas joining the variable selection
  beta_vs <- beta_value[group_info != 0]
  p <- length(beta_vs)
  eta <- 1/beta_tausq_value
  U <- runif(n = 1,
             min = 0,
             max = 1/(1+eta))
  ga.shape <- 0.5*p+0.5
  ga.rate <- 0.0
  num_group <- max(group_info)
  for (i in 1:num_group) {
    index <- group_info == i
    betai <- beta_value[index]
    ga.rate <- ga.rate + sum(betai^2/(2*beta_lambdasq_value[i]))
  }
  eta <- rtga(n = 1,
              shape = ga.shape,
              rate = ga.rate,
              a = 1,  # output will be between 0 and 1
              b = (1-U)/U)
  output <- 1/eta
  return(output)
}

update_nu <- function(y,
                      X,
                      intercept_value,
                      beta_a_value,
                      beta_b_value,
                      g_value,
                      dsq_value,
                      sigmasq_value,
                      delta_value,
                      nu_value,
                      n_core) {
  N <- length(y)
  logpdf <- function(log_nu_value) {
    output <- mst_reg_lpdf_parallel(y_vec_list = y,
                                    intercept = intercept_value,
                                    beta_a = beta_a_value,
                                    beta_b = beta_b_value,
                                    g_value_list = g_value,
                                    dsq = dsq_value,
                                    sigmasq = sigmasq_value,
                                    delta = delta_value,
                                    nu = exp(log_nu_value) + 2.0,  # the lower bound of nu is 2
                                    n_core = n_core)
    return(output)
  }
  
  output <- ESS(loglikelihood = logpdf,
                f.current = log(nu_value),
                mu = 0.0,
                Sigma = sqrt(1)^2)
  output <- exp(output) + 2.0
  return(output)
}

update_gof <- function(y,
                       g_value,
                       intercept_value,
                       beta_a_value,
                       beta_b_value,
                       B_value,
                       U_value,
                       delta_value,
                       dsq_value,
                       sigmasq_value,
                       nu_value,
                       b_scalar_value,
                       gof_K,
                       gof_L,
                       n_core) {
  N <- length(y)
  
  log_lik <- numeric(N)
  
  # step 1: define standardized residuals (wi)
  std.residual.list <- vector(mode = "list", length = N)
  loc.list <- vector(mode = "list", length = N)
  
  output.temp <- mclapply(1:N,
                          FUN = function(i) {
                            yi <- y[[i]]
                            gi <- g_value[[i]]
                            bi <- B_value[i]
                            ui <- U_value[i]
                            loc1i <- intercept_value + gi
                            loc2i <- beta_a_value + beta_b_value*loc1i
                            loci <- c(loc1i,loc2i)
                            residuali <- yi - loci - bi
                            std.residuali <- residuali/sqrt(sigmasq_value/ui)
                            std.residual.list <- std.residuali
                            loc.list <- loci + bi
                            log_lik <- mst_lpdf_c(x_vec = yi,
                                                  loc_value = loci,
                                                  dsq = dsq_value,
                                                  sigmasq = sigmasq_value,
                                                  delta = delta_value,
                                                  nu = nu_value)
                            output <- list(std.residual.list = std.residual.list,
                                           loc.list = loc.list,
                                           log_lik = log_lik)
                          },
                          mc.cores = n_core)
  
  for (i in 1:N) {
    log_lik[i] <- output.temp[[i]]$log_lik
    loc.list[[i]] <- output.temp[[i]]$loc.list
    std.residual.list[[i]] <- output.temp[[i]]$std.residual.list
  }
  
  std.residuals <- as.numeric(unlist(std.residual.list))
  
  location <- as.numeric(unlist(loc.list))
  
  # step 2: partition wi according to the values of location parameters
  thresholds <- qnorm(seq(0,1,length.out = gof_K + 1))
  
  # Use cut() to assign each w_i to a group based on locations values
  # 'groups' contains the group assignments for each w_i
  groups <- cut(location, breaks = thresholds, labels = FALSE, right = TRUE)
  
  # step 3: define chi square statistics dK
  # define the output of step 3 dK
  dK <- numeric(gof_K)
  
  # define the bin interval
  bin.interval <- seq(0-1e-5,1,length.out = gof_L + 1)
  
  # calculate dK for k = 1,...,K
  for (k in 1:gof_K) {
    index <- groups == k
    if (sum(index)) {
      # find the standardized residuals in group k
      std.residuals.k <- std.residuals[index]
      
      # the quantity used in the cut function
      quantity.k <- pnorm(std.residuals.k)
      
      # determine the cell for group k
      cell.k <- cut(quantity.k, breaks = bin.interval, labels = FALSE, right = TRUE)
      
      # nk: the total number of residuals in group k
      nk <- length(std.residuals.k)
      dkL <- numeric(gof_L)
      # calculate dkl
      for (l in 1:gof_L) {
        # Ol: the observed number of residuals in bin l
        Ol <- sum(cell.k == l)
        # pl <- pnorm(bin.interval[l+1]) - pnorm(bin.interval[l])
        pl <- bin.interval[l+1] - bin.interval[l]
        dkL[l] <- ((Ol - nk*pl)/sqrt(nk*pl))^2
      }
      dK[k] <- sum(dkL) 
    }
  }
  # step 4: compute the test statistic d
  d <- sum(dK)
  
  output <- list(d = d,
                 log_lik = log_lik)
  
  return(output)
}
