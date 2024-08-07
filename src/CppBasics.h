using namespace Rcpp;

// Copy elements of one vector containing real numbers to another vector
// from : vector whose elements are to be copied
// to : vector where the copied elements are stored
void vec_copy(NumericVector from, NumericVector to){
  for(int i=0; i<to.size(); i++)
    to[i] = from[i] ;
}

// Calculate the p-norm of a vector as ||x||_p = (\sum |x_i|^p)^(1/p)
// x : vector whose norm is to be calculated
// p : integer denoting the p-th norm
// power : boolean; If false, ||x||_p is returned, else (||x||_p)^p is returned
double norm_c(NumericVector x, int p=2, bool power=false){
  double res=0.0 ;

  // res = \sum |x_i|^p
  for(int i=0; i<x.size(); i++)
    res += pow(abs(x[i]), p) ;
  if(power == true)
    return(res) ;   //if power = true, return \sum |x_i|^p
  else
    return(pow(res, 1.0/(double)p)) ;   // otherwise return (\sum |x_i|^p)^(1/p)
}

// Multiply a matrix with a vector and store it in "res" vector passed as an argument
// A : matrix
// x : vector
// res : vector where A %*% x is stored
void product_vec_c(NumericMatrix A, NumericVector x, NumericVector res){
for(int i=0; i<A.nrow(); i++)	{
  res[i] = 0.0 ;
  for(int j=0; j<A.ncol(); j++)
    res[i] += A(i,j)*x[j] ; // res = A*x
  }
}

// Multiply a matrix with a scalar and store it in "res" matrix passed as an argument
// A : matrix
// x : scalar
// res : matrix where A * x is stored
void product_c(NumericMatrix A, double x, NumericMatrix res){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      res(i,j) = x * A(i,j) ;   // Multiply every element of A with x
}

// Calculate the transpose of a matrix and store in a matrix passed as an argument
// A : matrix whose transpose is to be calculated
// t_A : matrix where t(A) is stored
void transpose_c(NumericMatrix A, NumericMatrix t_A){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      t_A(j,i) = A(i,j) ;   // calculate transpose of A, i.e. A_ij gets stored as t(A)_ji
}

// Calculate Cholesky Decomposition of a matrix and store in a matrix passed as an argument
// A : matrix whose cholesky decomposition is required
// res : matrix where the decomposition is stored
void chol_c(NumericMatrix A, NumericMatrix res){
  int n = A.nrow() ;
  arma::mat arma_A;
  arma_A = as<arma::mat>(A);
  arma::mat tmp;
  if(arma_A.is_symmetric()){
    tmp = arma::chol(arma_A) ;    // Cholesky decomposition using Armadillo
  }
  else{
    arma_A = arma::symmatu(arma_A); // make arma_A a symmetric matrix
    tmp = arma::chol(arma_A) ;
  }
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      res(i,j) = tmp(i,j) ;     // storing the result in res
}

// Calculate the inverse of a matrix and store in a matrix passed as an argument
// A : matrix whose inverse is to be calculated
// A_inv : matrix where inverse(A) is stored
void solve_store_c(NumericMatrix A, NumericMatrix A_inv){
  arma::mat tmp = as<arma::mat>(A).i() ; // Calculate inverse(A) using Armadillo

  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      A_inv(i,j) = tmp(i,j) ;     // store the inverse in res
}

// Calculate the inverse of a matrix and return it
// A : matrix whose inverse is to be calculated
NumericMatrix solve_c(NumericMatrix A){
  arma::mat res = as<arma::mat>(A).i() ;    // Calculate inverse(A) using Armadillo
  return(wrap(res)) ;
}

// Calculate A^{\top} * x, A: matrix and X: vector
// A : matrix
// x : vector
// res : vector where transpose(A) %*% x is stored
void Atx_c(NumericMatrix A, NumericVector x, NumericVector res){
  for(int i=0; i<A.ncol(); i++)	{
    res[i] = 0.0 ;
    for(int j=0; j<A.nrow(); j++)
      res[i] += A(j,i)*x[j] ;   // Calculate res = transpose(A)*x
  }
}

// Calculate A^{\top} * A, where A is a matrix
// A : matrix
// res : matrix where transpose(A) %*% A is stored
void AtA_c(NumericMatrix A, NumericMatrix res){
  // res = t(A) %*% A
  for(int i=0; i<A.ncol(); i++){
    for(int j=0; j<A.ncol(); j++){
      res(i,j) = 0.0 ;
      for(int k=0; k<A.nrow(); k++)
        res(i,j) += A(k,i)*A(k,j) ;   // Calculating transpose(A) %*% A element wise
    }
  }
}

// Calculate the innerproduct of 2 vectors x and y
// x, y : vectors whose inner product need to be calculated
double innerProduct_c(NumericVector x, NumericVector y){
  double res=0.0 ;
  for(int i=0; i<x.size(); i++)
    res += x[i]*y[i] ;    // res stores x*y
  return(res) ;
}

// Draw a sample from multivariate Normal (mu, Sigma) and return it
// mu : mean vector of the Normal Distribution
// Sigma : variance vector of the Normal Distribution
// res : vector that stores the drawn sample from multivariate Normal (mu, Sigma)
void rmvnorm_c(NumericVector mu, NumericMatrix Sigma, NumericVector res){
  // Use Armadillo to generate a sample from Normal (0, I). Call it Y.
  arma::mat Y = arma::randn(1,mu.size());

  // Z  = \mu + (\Sigma^(1/2)) * Y is a sample from Normal (mu, Sigma)
  arma::mat Z = arma::repmat(as<arma::vec>(mu),1,1).t() + Y * arma::chol(as<arma::mat>(Sigma)) ;
  for(int i=0; i<res.size(); i++)
    res[i] = Z(0,i) ;
}

// Add 2 matrices and store it in a matrix passed as an argument
// A, B : matrices to be added
// res : matrix that stores the sum of A and B
void sum_c(NumericMatrix A, NumericMatrix B, NumericMatrix res){
  for(int i=0; i<A.nrow(); i++)
    for(int j=0; j<A.ncol(); j++)
      res(i,j) = A(i,j) + B(i,j) ;    // Add A and B element wise
}
