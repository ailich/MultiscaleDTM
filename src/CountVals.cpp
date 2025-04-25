#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector C_CountVals(NumericVector z, size_t ni, size_t nw) {
  NumericVector out(ni, NA_REAL);
  const double* z_ptr = z.begin(); // Direct pointer access to z
  
  for (size_t i = 0; i < ni; ++i) {
    size_t start = i * nw;
    size_t count = 0;
    for (size_t j = 0; j < nw; ++j) {
      double val = z_ptr[start + j];
      if (!std::isnan(val)) {
        ++count;
      }
    }
    
    out[i] = count;
  }
  
  return out;
}
