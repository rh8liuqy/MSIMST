#include <RcppArmadillo.h>
#include <math.h>
#include "CppBasics.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// The function to calculate one component of the phiX matrix.
double psi(double x){
  if( x < -1.0 )
    return(0.0) ;   // return 0 for x < -1
  else if( x < 0.0 )
    return( 0.5*(x+1.0)*(x+1.0) ) ;   // return 1/2*(x+1)^2 if -1<x<0
  else if(x < 1.0)
    return( 0.5 + 0.5*(2.0-x)*x ) ;   // return 1/2 * 1/2*(2-x)x if 0<=x<1
  else
    return(1.0) ;   // return 1 if x >= 1
}

// psi_l(x) = \int_{-1}^x h_l(t) dt
// l : index ranging from 0 to (#knots - 1)
// u : specified knots from -1 to 1
double psi_l(double x, int l, const NumericVector& u){
  if(l==0){
    // specific algerba for psi_0
    if(x <= -1)
      return (0.0);
    else if(x < u[1]){
      return( ( u[1]*(1+x) + 0.5 - 0.5*x*x )/(u[1]+1) ) ;
    }
    else
      return(0.5*(u[1]+1)) ;
  }
  else{
    double u1, tmp ;
    u1 = u[l-1] ;
    
    if(x < u[l]){
      // if x < u_l, return (u_l - u_{l-1})* \psi( (x - u_l)/(u_l - u_{l-1}) )
      
      tmp = u[l] - u1 ;
      return( tmp*psi((x-u[l])/tmp) ) ;
    }
    else{
      // if x >= u_l, return 0.5(u_l - u_{l-1})* (u_{l-1} - u_l){\psi( (x - u_l)/(u_{l-1} - u_l) ) - 0.5}
      
      double u3 ;
      int L=u.size()-1 ;
      if(l==L)
        u3 = 2.0*u[L] - u[L-1] ;
      else
        u3 = u[l+1] ;
      tmp = u3 - u[l] ;
      return( 0.5*(u[l]-u1) + tmp*(psi((x-u[l])/tmp)-0.5) ) ;
    }
  }
}

//' The Function to Calculate the phiX Matrix for Estimating Single-Index Function
//' 
//' @param Xbeta The single index values. A vector of length n.
//' @param u The vector spanning from -1 to 1 with length L + 1.
//' @param L An integer defining the number of nodes.
//' @return A n by L + 1 matrix.
//' @examples
//' L <- 50
//' u <- seq(-1,1,length.out = L + 1)
//' phiX <- phiX_c(0.5,u,L)
//' print(phiX)
//' @description
//' The function \verb{phiX_c} is used to generate the phiX matrix associated with the Gaussian process prior.
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix phiX_c(const NumericVector& Xbeta,const NumericVector& u,int L){
  int n = Xbeta.size();
  NumericMatrix Psi(n,L+1);
  for(int i=0; i<n; i++){
    for(int l=0; l<L+1; l++){
      Psi(i,l) = psi_l(Xbeta[i], l, u) ;
    }
  }
  return(Psi);
}

// Function to calculate cs_inv
// [[Rcpp::export]]
arma::mat cs_inv_c(double dsq, double sigmasq, int n) {
  arma::mat p1 = arma::eye<arma::mat>(n, n) / sigmasq;
  double common_term = dsq / (sigmasq * sigmasq + sigmasq * dsq * n);
  arma::mat p2(n, n, arma::fill::none);
  p2.fill(common_term);
  
  arma::mat output = p1 - p2;
  return output;
}

// Function to calculate b.func
// [[Rcpp::export]]
double b_func_c(double nu) {
  double output = -sqrt(nu / M_PI) * R::gammafn(0.5 * (nu - 1)) / R::gammafn(0.5 * nu);
  return output;
}

// Function to calculate cs_det
// [[Rcpp::export]]
double cs_det_c(double dsq, double sigmasq, int n) {
  double p1 = n * dsq * pow(sigmasq, n - 1);
  double p2 = pow(sigmasq, n);
  double output = p1 + p2;
  return output;
}

// Function to calculate the log pdf of multivariate Student t distribution for one observation
// [[Rcpp::export]]
double mvt_lpdf_c(NumericVector x_vec, NumericVector mu_vec, arma::mat Sigma_inv, double Sigma_det, double nu) {
  int p = x_vec.size();
  double p1 = lgamma(0.5 * (nu + p));
  double p2 = lgamma(0.5 * nu) + 0.5 * p * log(nu) + 0.5 * p * log(M_PI) + 0.5 * log(Sigma_det);
  
  arma::vec xstar = as<arma::vec>(x_vec) - as<arma::vec>(mu_vec);
  double dx = as_scalar(xstar.t() * Sigma_inv * xstar);
  double p3 = -0.5 * (nu + p) * log(1 + (1/nu) * dx);
  
  double output = p1 - p2 + p3;
  return output;
}

// Function to calculate mst_lpdf
// [[Rcpp::export]]
double mst_lpdf_c(NumericVector x_vec, NumericVector loc_value, 
                  double dsq, double sigmasq, double delta, double nu) {
  int p = x_vec.size();
  double b_scalar = b_func_c(nu);
  NumericVector mu_vec = loc_value + b_scalar * delta;
  double p1 = log(2);
  // lpdf of the multivariate student t 
  arma::mat Sigma_inv = cs_inv_c(dsq + delta * delta, sigmasq, p);
  double Sigma_det = cs_det_c(dsq + delta * delta, sigmasq, p);
  double p2 = mvt_lpdf_c(x_vec, mu_vec, Sigma_inv, Sigma_det, nu);
  // lcdf of the univariate student
  double lambda = 1 - delta * delta * accu(Sigma_inv);
  arma::vec x_star = x_vec - mu_vec;
  arma::vec Sigma_inv_x_star = Sigma_inv * x_star;
  double dx = arma::as_scalar(x_star.t() * Sigma_inv_x_star);
  double lcdf_ub = delta * accu(Sigma_inv_x_star) * sqrt((nu + p) / (nu + dx));
  lcdf_ub = lcdf_ub / sqrt(lambda);
  double p3 = R::pt(lcdf_ub, nu + p, true, true);
  double output = p1 + p2 + p3;
  return output;
}

// calculates B_{L,j}(t)
// L - degree of bernstein polynomials, number of basis functions: L+1
// j - indicates which of the L+1 bernstein polynomials are being evaluated
// t - the value the polynomial is being evaluated at, in [-1,1]
double Btilde(int L, int j, double t) {
  
  // if(t >= 1.0 || t <= -1.0)
  //   return(0.0) ;
  // else
  return((0.5/(L+1)) * R::dbeta((0.5*(t+1)), (j+1), (L-j+1), false)) ;
}

// calculates D_beta matrix
// t is a vector of length n, t in [-1,1]
// L is the degree of the bernstein polynomials
//[[Rcpp::export]]
NumericMatrix Dbeta(int L, NumericVector t) {
  
  int n = t.size();
  NumericMatrix BigB(n, L), A(L, L), D_mat(n, L);
  
  for(int i=0; i<n; i++)
    for(int j=0; j<L; j++)
      BigB(i,j) = Btilde(L, j+1, t[i]) ;
  // rows are over t, columns are over 0:L, values[i,j] is Btilde(L,j,t[i])
  
  for(int i=0; i<L; i++)
    for(int j=0; j<L; j++){
      if(i >= j) {
        A(i, j) = 1;
      }
      else {
        A(i, j) = 0;
      }
    }
    
  arma::mat BigB_mat = as<arma::mat>(BigB);
  arma::mat A_mat = as<arma::mat>(A);
  arma::mat D_armamat = as<arma::mat>(D_mat);
  D_armamat = BigB_mat * A_mat;
  D_mat = as<NumericMatrix>(wrap(D_armamat));
  return(D_mat);
}

// Function to draw a sample of size n from truncated Normal distribution using an exact Hamiltonian Monte Carlo Algorithm
// Thank Snigdha Das from Texas A&M University for sharing codes of her implementation.
// x_init ~ d-variate Normal(mu, Sigma)
// constraints: <ff(j,_), x_init> + gg[j] \geq 0, j=1,...,m
// n: required sample size
// mu, Sigma : Parameters of the normal distribution
// x_init : initial value of the observation that follows the d-variate Normal distribution
// ff, gg : functions that define constraints
// n_burn : Burn in period of the Monte Carlo Algorithm
// [[Rcpp::export]]
NumericMatrix rtmvnormHMC(int n, const NumericVector& mu, const NumericMatrix& Sigma,
                          const NumericVector& x_init, const NumericMatrix& ff,
                          const NumericVector& gg, int n_burn){
  
  // Initialize required varaibles
  int d=mu.size(), m=ff.nrow(), h ;
  double u, phi, tmp, tmp2, T_end, T_h, alpha_h ;
  NumericVector x(d), s(d), a(d), b(d), T(m), x_dot(d), g(m) ;
  NumericMatrix f(m,d), res(n,d), Sigma_chol_L(d,d), Sigma_chol_U(d,d) ;
  
  // Sigma = Sigma_col_L * Sigma_col_U by Cholesky Decomposition
  chol_c(Sigma, Sigma_chol_U) ;
  transpose_c(Sigma_chol_U, Sigma_chol_L) ;
  
  // x = (Sigma_chol_L)^-1 (x_init-mu) ~ d-variate Normal(0, I_d)
  product_vec_c(solve_c(Sigma_chol_L), x_init-mu, x) ;
  
  // Adjust the constraints on x_init to apply them on x
  for(int j=0; j<m; j++){
    g[j] = innerProduct_c(ff(j,_), mu) + gg[j] ;
    for(int i=0; i<d; i++)
      f(j,i) = innerProduct_c(Sigma_chol_U(i,_), ff(j,_)) ;
  }
  
  // Perform the HMC algorithm
  for(int nn=0; nn<n+n_burn; nn++){
    
    s = rnorm(d) ;    // Draw s ~ Normal(0, I_d)
    vec_copy(s, a) ;   // a = s
    vec_copy(x, b) ;   // b = s
    
    T_end = M_PI/2.0 ;    // T_end = PI/2, end time point
    
    while(1){
      for(int j=0; j<m; j++){
        
        // u = sqrt((<f_j,a>)^2 + (<f_j,b>)^2)
        tmp = innerProduct_c(f(j,_), a) ;
        tmp2 = innerProduct_c(f(j,_), b) ;
        u = sqrt(tmp*tmp + tmp2*tmp2) ;
        
        // phi = tan^-1 (-tmp/tmp2)
        phi = atan2(-1.0*tmp, tmp2) ;
        
        if((u < g[j]) || (u < -1.0*g[j]))     // condition to check if the jth particle has hit the wall
          T[j] = T_end ;    // Set the time point of the jth particle as the end point
        else{
          // if the jth particle is still along trajectory, calculate the time to hit the wall using
          // u * cos(T_j + phi) + g_j = 0
          T[j] = acos(-1.0*g[j]/u) - phi ;
          if(T[j] < 0.0){
            T[j] += 2.0*M_PI ;
            // Rprintf("1\n") ;
          }
        }
      }
      
      // T_h = min {T_1, T_2, ..., T_m}
      h = 0 ;
      T_h = T[0] ;
      for(int j=1; j<m; j++){
        if(T[j] < T_h){
          T_h = T[j] ;
          h = j ;
        }
      }
      
      if(T_h < T_end){
        // if atleast one of the particles has not hit the wall, Set:
        // x_i = a_i sin(T_h) + b_i cos(T_h)
        // x_dot_i = a_i sin(T_h) - b_i cos(T_h)
        for(int i=0; i<d; i++){
          tmp = sin(T_h - 1e-10) ;
          tmp2 = cos(T_h - 1e-10) ;
          x[i] = a[i]*tmp + b[i]*tmp2 ;
          x_dot[i] = a[i]*tmp2 - b[i]*tmp ;
        }
        
        // set a_i as the reflected velocity at time T_h
        alpha_h = innerProduct_c(f(h,_), x_dot) / norm_c(f(h,_),2,true) ;
        for(int i=0; i<d; i++)
          a[i] = x_dot[i] - 2.0*alpha_h*f(h,i) ;
        
        // Set b_i = x_i
        vec_copy(x, b) ;
        
        // Decrease the end time point by T_h
        T_end -= T_h ;
      }
      else{
        
        // if all particles have hit the wall, Set x_i = a_i sin(T_h) + b_i cos(T_h)
        // and stop the algorithm
        for(int i=0; i<d; i++){
          tmp = sin(T_end) ;
          tmp2 = cos(T_end) ;
          x[i] = a[i]*tmp + b[i]*tmp2 ;
        }
        break ;
      }
    }
    // Store samples after burn-in period
    if(nn >= n_burn)
      for(int i=0; i<d; i++)
        res(nn-n_burn,i) = innerProduct_c(Sigma_chol_L(i,_), x) + mu[i] ;
    // Obtained x is a sample from truncated Normal(0, I_d)
    // Transform it to get the required truncated Normal(mu, Sigma)
  }
  return(res) ;
}
