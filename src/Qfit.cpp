#define ARMA_WARN_LEVEL 1
#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//Fit Wood/Evans Quadratic Surface or Plane with Intercept

//na.rm=TRUE, force_center=FALSE
// [[Rcpp::export]]
arma::mat C_Qfit1_narmT(const arma::vec& z,
                             const arma::mat& X_full,
                             const bool return_intercept,
                             size_t ni, size_t nw) {
  
  unsigned int thresh = X_full.n_cols;
  
  int nlyr = return_intercept ? X_full.n_cols : X_full.n_cols - 1;
  
  arma::mat out(ni, nlyr, arma::fill::value(NA_REAL));
  
  for (size_t i = 0; i < ni; ++i) {
    size_t start = i * nw;
    arma::vec zw = z.subvec(start, start + nw - 1);
    arma::uvec valid_idx = arma::find_finite(zw);
    
    if (valid_idx.n_elem >= thresh) {
      arma::vec zw_clean = zw.elem(valid_idx);
      arma::mat X = X_full.rows(valid_idx);
      
      arma::vec unique_vals = arma::unique(zw_clean);
      if (unique_vals.n_elem == 1) {
        out.row(i).zeros();
        if (return_intercept) {
          out(i, nlyr - 1) = unique_vals[0];  // Last coefficient (e.g., intercept "f")
          }
        } else {
        // Solve using least squares: B = solve(X, Z)
        arma::vec coef = arma::solve(X, zw_clean);
          if (return_intercept) {
            out.row(i) = coef.t();  // Full coefficients (including intercept)
          } else {
            out.row(i) = coef.subvec(0, coef.n_elem - 2).t();  // Exclude intercept
          }}
    }
  }
  return out;
}

//na.rm=FALSE, force_center=FALSE
// [[Rcpp::export]]
arma::mat C_Qfit1_narmF(const arma::vec& z,
                             const arma::mat& X,
                             const arma::mat& Xt,
                             const arma::mat& XtX_inv,
                             const bool return_intercept,
                             size_t ni, size_t nw) {
  
  int nlyr = return_intercept ? X.n_cols : X.n_cols - 1;
  arma::mat out(ni, nlyr, arma::fill::value(NA_REAL));
  
  for (size_t i = 0; i < ni; ++i) {
    size_t start = i * nw;
    arma::vec zw = z.subvec(start, start + nw - 1);
    
    if (!zw.has_nan()) {
      arma::vec unique_vals = arma::unique(zw);
      
      if (unique_vals.n_elem == 1) {
        out.row(i).zeros();
        if (return_intercept) {
          out(i, nlyr - 1) = unique_vals[0];  // Last coefficient (e.g., intercept "f")
        }
      } else {
        arma::vec coef = XtX_inv * (Xt * zw);
          if (return_intercept) {
            out.row(i) = coef.t();  // Full coefficients (including intercept)
          } else {
            out.row(i) = coef.subvec(0, coef.n_elem - 2).t();  // Exclude intercept
          }}
    }
  }
  
  return out;
}


//Fit Wood/Evans Quadratic Surface or Plane forced through center

//na.rm=TRUE, force_center=TRUE
// [[Rcpp::export]]
arma::mat C_Qfit2_narmT(const arma::vec& z,
                   const arma::mat& X_full,
                   size_t ni, size_t nw) {
  
  int nlyr = X_full.n_cols;
  arma::mat out(ni, nlyr, arma::fill::value(NA_REAL));
  
  unsigned int thresh = nlyr; //Setting to nlayer makes it work for planar fit too
  
  for (size_t i = 0; i < ni; ++i) {
    size_t start = i * nw;
    arma::vec zw = z.subvec(start, start + nw - 1);
    
    // Centering the window values
    double center_val = zw(nw / 2);
    zw -= center_val;
    
    arma::uvec valid_idx = arma::find_finite(zw);
    if (valid_idx.n_elem >= thresh) {
      arma::vec zw_clean = zw.elem(valid_idx);
      arma::mat X = X_full.rows(valid_idx);
      
      arma::vec unique_vals = arma::unique(zw_clean);
      if (unique_vals.n_elem == 1) {
        out.row(i).zeros();
        } else {
        // Least squares solution
        arma::vec coef = arma::solve(X, zw_clean);
        out.row(i) = coef.t(); // row vector
      }
      }
    }
  
  return out;
}

//na.rm=FALSE, force_center=TRUE
// [[Rcpp::export]]
arma::mat C_Qfit2_narmF(const arma::vec& z,
                             const arma::mat& X,
                             const arma::mat& Xt,
                             const arma::mat& XtX_inv,
                             size_t ni, size_t nw) {
  
  int nlyr = X.n_cols;
  arma::mat out(ni, nlyr, arma::fill::value(NA_REAL));
  
  for (size_t i = 0; i < ni; ++i) {
    size_t start = i * nw;
    arma::vec Z = z.subvec(start, start + nw - 1);
    
    // Center Z using middle value
    double center_val = Z(Z.n_elem / 2);
    Z -= center_val;
    
    if (!Z.has_nan()) {
      arma::vec unique_vals = arma::unique(Z);
      
      if (unique_vals.n_elem == 1) {
        out.row(i).zeros();  // All coefficients = 0
      } else {
        // OLS coefficients using precomputed XtX_inv and Xt
        arma::vec coef = XtX_inv * (Xt * Z);
        out.row(i) = coef.t();
      }
    }
    }
  
  return out;
}


//Planar fit (Not needed. Generalized Qfit)

//na.rm=TRUE, force_center=FALSE
// arma::mat C_Pfit1_narmT(const arma::vec& z,
//                         const arma::mat& X_full,
//                         size_t ni, size_t nw) {
//   
//   int nlyr = X_full.n_cols-1;
//   arma::mat out(ni, nlyr, arma::fill::value(NA_REAL));
//   
//   unsigned int thresh = 3;
//   
//   for (size_t i = 0; i < ni; ++i) {
//     size_t start = i * nw;
//     arma::vec zw = z.subvec(start, start + nw - 1);
//     arma::uvec valid_idx = arma::find_finite(zw);
//     
//     if (valid_idx.n_elem >= thresh) {
//       arma::vec zw_clean = zw.elem(valid_idx);
//       arma::mat X = X_full.rows(valid_idx);
//       
//       arma::vec unique_vals = arma::unique(zw_clean);
//       if (unique_vals.n_elem == 1) {
//         out.row(i).zeros();
//       } else {
//         // Solve using least squares: B = solve(X, Z)
//         arma::vec coef = arma::solve(X, zw_clean);
//         out.row(i) = coef.subvec(0, coef.n_elem - 2).t();  // Transpose to write as row
//       }
//     }
//   }
//   
//   return out;
// }

//na.rm=FALSE, force_center=FALSE (Not Used)
// arma::mat C_Pfit1_narmF(const arma::vec& z,
//                         const arma::mat& X,
//                         const arma::mat& Xt,
//                         const arma::mat& XtX_inv,
//                         size_t ni, size_t nw) {
//   
//   int nlyr = X.n_cols - 1;
//   arma::mat out(ni, nlyr, arma::fill::value(NA_REAL));
//   
//   for (size_t i = 0; i < ni; ++i) {
//     size_t start = i * nw;
//     arma::vec zw = z.subvec(start, start + nw - 1);
//     
//     if (!zw.has_nan()) {
//       arma::vec unique_vals = arma::unique(zw);
//       
//       if (unique_vals.n_elem == 1) {
//         out.row(i).zeros();
//       } else {
//         arma::vec coef = XtX_inv * (Xt * zw);
//         out.row(i) = coef.subvec(0, coef.n_elem - 2).t();  // Transpose to write as row
//       }
//     }
//   }
//   return out;
// }

