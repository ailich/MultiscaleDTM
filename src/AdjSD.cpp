#define ARMA_WARN_LEVEL 1
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//SD from planar fit with na.rm=TRUE
// [[Rcpp::export]]
arma::vec C_AdjSD_narmT(const arma::vec& z, 
                        const arma::mat& X_full, 
                        size_t ni, size_t nw) {
  arma::vec out(ni, arma::fill::value(NA_REAL));
  unsigned int thresh = 4;
  
  for (size_t i = 0; i < ni; ++i) {
    size_t start = i * nw;
    arma::vec zw_full = z.subvec(start, start + nw - 1);
    arma::uvec non_na_idx = arma::find_finite(zw_full);
    
    if (non_na_idx.n_elem >= thresh) {
      arma::vec zw = zw_full.elem(non_na_idx);
      
      arma::vec unique_vals = arma::unique(zw);
      if (unique_vals.n_elem == 1) {
        out[i] = 0;
      } else {
        arma::mat X = X_full.rows(non_na_idx);
        arma::vec resid = zw - X * arma::solve(X, zw);
        out[i] = arma::stddev(resid);
      }
    }
  }
  return out;
}

//SD from planar fit with na.rm=FALSE
// [[Rcpp::export]]
arma::vec C_AdjSD_narmF(const arma::vec& z, 
                        const arma::mat& X, 
                        const arma::mat& Xt, 
                        const arma::mat& XtX_inv, 
                        size_t ni, size_t nw) {
  arma::vec out(ni, arma::fill::value(NA_REAL));
  
  for (size_t i = 0; i < ni; ++i) {
    size_t start = i * nw;
    arma::vec Z = z.subvec(start, start + nw - 1);
    
    if (!Z.has_nan()) {
      arma::vec unique_vals = arma::unique(Z);
      if (unique_vals.n_elem == 1) {
        out[i] = 0;
      } else {
        arma::vec resid = Z - X * (XtX_inv * (Xt * Z));
        out[i] = arma::stddev(resid);
      }
    }
  }
  
  return out;
}