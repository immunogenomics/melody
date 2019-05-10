#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat moe_gene(const arma::mat & X, const arma::mat & Phi, const arma::mat & R, float lambda, const uvec & idx_numerical) {
  arma::mat X_corr = X;
  unsigned K = R.n_rows;
    
  arma::mat Phi_Rk, W;
  arma::vec lambda_vec(Phi.n_rows);
  lambda_vec.fill(lambda);
  lambda_vec.row(0).zeros(); // no regularization on intercept
//   lambda_vec.rows(idx_numerical).zeros(); // no regularization on intercept
  arma::mat lambda_mat = arma::diagmat(lambda_vec);

  arma::mat mu = arma::zeros<arma::mat>(X.n_rows, X.n_cols);

  // B (numerical values only) x K
  arma::mat A = Phi.rows(idx_numerical) * R.t() * arma::diagmat(1 / arma::sum(R, 1));
  for (unsigned k = 0; k < K; k++) {
      
    // PT1: compute model parameters
    // this is a log-linear model: log(y) ~ beta
    Phi_Rk = Phi * arma::diagmat(R.row(k));
    W = arma::inv(Phi_Rk * Phi.t() + lambda_mat) * Phi_Rk * arma::log1p(X.t()); 
    
    // PT2: compute cell residuals
    // residuals should be additive: y = exp(beta) + resid
//     W.row(0).zeros(); 
    mu = W.t() * Phi;
    X_corr -= expm1(mu);
    
    // PT3: add back mean value for continuous value
    // use cols and rows so A and R stay matrices
    X_corr += expm1(W.rows(idx_numerical).t() * A.cols(k, k) * R.rows(k, k));
  }

  // PT3: add residuals to model means
  
  return X_corr;
}
// [[Rcpp::export]]
Rcpp::List foo(const arma::mat & X, const arma::mat & Phi, const arma::mat & R, float lambda, const uvec & idx_numerical) {
  unsigned K = R.n_rows;
    
  arma::mat Phi_Rk;
  Rcpp::List W(K);
  arma::vec lambda_vec(Phi.n_rows);
  lambda_vec.fill(lambda);
  lambda_vec.rows(idx_numerical).zeros(); // no regularization on intercept
  arma::mat lambda_mat = arma::diagmat(lambda_vec);

  arma::mat mu = arma::zeros<arma::mat>(X.n_rows, X.n_cols);

  // B (numerical values only) x K
  arma::mat A = Phi.rows(idx_numerical) * R.t() * arma::diagmat(1 / arma::sum(R, 1));
  for (unsigned k = 0; k < K; k++) {
      
    // PT1: compute model parameters
    // this is a log-linear model: log(y) ~ beta
    Phi_Rk = Phi * arma::diagmat(R.row(k));
    W(k) = arma::inv(Phi_Rk * Phi.t() + lambda_mat) * Phi_Rk * arma::log1p(X.t()); 
    
  }  
  return W;
}

// // [[Rcpp::export]]
// Rcpp::List get_model_pars_rmoe(const arma::mat & X, const arma::mat & Phi, const arma::mat & R, float lambda) {    
//   arma::mat Phi_Rk;
//   arma::vec lambda_vec(Phi.n_rows);
//   lambda_vec.fill(lambda);
//   arma::mat lambda_mat = arma::diagmat(lambda_vec);
    
//   Rcpp::List W(R.n_rows);
//   for (int k = 0; k < R.n_rows; k++) {
//     Phi_Rk = Phi * arma::diagmat(R.row(k));
//     W[k] = arma::inv(Phi_Rk * Phi.t() + lambda_mat) * Phi_Rk * X.t();
//   }
//   return W;
// }
