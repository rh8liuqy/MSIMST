#' Weighted Finite Population Bayesian Bootstrap
#'
#' @param y The index of survey data.
#' @param w Survey weights. The summation of survey weights should equal the population size
#' @param N The population size.
#' @param n The sample size.
#'
#' @return The re-sampled index of survey data.
#' @description
#' The function is implemented based on the WFPBB algorithm from \insertCite{Gunawan2020}{MSIMST}.
#' 
#' @references 
#' \insertRef{Gunawan2020}{MSIMST}
#' @examples
#' set.seed(100)
#' output_data <- reg_simulation3(N = 5000,
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
#' # set the population size
#' population_N <- 5000
#' # set the sample size
#' n <- length(y)
#' # run the WFPBB algorithm
#' index_WFPBB <- WFPBB(y = 1:n,
#'                      w = survey_weight,
#'                      N = population_N,
#'                      n = n)
#' print(head(index_WFPBB))
#' 
#' @export
#' 
WFPBB <- function(y,
                  w,
                  N,
                  n) {
  if (length(y) != n) {
    stop("length(y) must equal n !")
  }
  if (length(w) != n) {
    stop("length(w) must equal n !")
  }
  
  # l: number of bootstrap selections
  l <- numeric(n)
  
  Nstar <- (N-n)/n
  
  # output: the output of the for loop
  output <- numeric(N - n)
  for (k in 1:(N-n)) {
    # print the log information if the population size is large
    if (k %% 10000 == 1) {
      print(paste0("WFPBB procedure || iteration: ", k, " of ", N - n))
    }
    probk <- (w-1+l*Nstar)/(N-n+(k-1)*Nstar)
    ykstar <- sample(x = y, size = 1, replace = FALSE, prob = probk)
    output[k] <- ykstar
    l[ykstar] <- l[ykstar] + 1
  }
  
  # stacking to form a pseudo population
  output <- c(y,output)
  
  # randomly draw a sample of size n from the pseudo population
  output <- sample(x = output, size = n, replace = FALSE)
  return(output)
}
