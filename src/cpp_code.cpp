#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//Extracts relevant window from matrix based on position of central pixel and window size
// [[Rcpp::export]]
NumericMatrix C_extract_window(NumericMatrix r, IntegerVector w, IntegerVector idx){
  int nr = w(0);
  int nc = w(1);
  
  int rast_row_center = idx(0);
  int rast_row_top = rast_row_center-(nr-1)/2;
  
  int rast_col_center= idx(1);
  int rast_col_left = rast_col_center-(nc-1)/2;
  
  IntegerMatrix r_idx(nr, nc);
  IntegerMatrix c_idx(nr, nc);
  for(int i=0; i < nr; ++i){
    for(int j=0; j < nc; ++j){
      r_idx(i,j) = i+rast_row_top;
      c_idx(i,j)=j+rast_col_left;
    }}
  NumericMatrix dat(nr, nc);
  for(int k=0; k < dat.size(); ++k){
    dat[k]= r(r_idx[k], c_idx[k]);
  }
  return dat; //extracted window as a matrix
}


//Subsets rows in a numeric mattrix according to a logical vector
// [[Rcpp::export]]
NumericMatrix subset_mat_rows(NumericMatrix x, LogicalVector idx) {
  NumericVector x2 = as<NumericVector>(x);
  LogicalVector idx2 = rep(idx, x.ncol());
  NumericVector out_vect = x2[idx2];
  NumericMatrix out_mat(sum(idx), x.ncol(), out_vect.begin());
  return out_mat;
}

//Ordinary Least Squares (only returns parameters)
// [[Rcpp::export]]
NumericVector C_OLS_params(arma::mat X, arma::mat Y){
  arma::mat Xt = trans(X);
  arma::mat XtX = Xt * X;
  double d = det(XtX);
  if(d==0){
    NumericVector B2(X.n_cols, NA_REAL);
    return B2;
  } else{
    arma::mat XtX_inv= inv(XtX);
    arma::mat H = X * XtX_inv * Xt;
    NumericVector B= Rcpp::as<Rcpp::NumericVector>(wrap(XtX_inv * (Xt * Y)));
    return B;
  }}

//Ordinary Least Squares (only returns residuals)
// [[Rcpp::export]]
NumericVector C_OLS_resid(arma::mat X, arma::mat Y){
  arma::mat Xt = trans(X);
  arma::mat XtX = Xt * X;
  double d = det(XtX);
  if(d==0){
    NumericVector resid2(X.n_rows, NA_REAL);
    return resid2;
  } else{
    arma::mat XtX_inv= inv(XtX);
    arma::mat H = X * XtX_inv * Xt;
    arma::mat Yhat = H * Y;
    NumericVector resid = Rcpp::as<Rcpp::NumericVector>(wrap(Yhat - Y));
    return resid;
    }
}

//Multiscale metrics across matrix using sliding window (Quadratic Fit)
// [[Rcpp::export]]
NumericMatrix C_multiscale2(NumericMatrix r, IntegerVector w, NumericMatrix X, bool na_rm){
  int nr= r.nrow();
  int nc= r.ncol();
  int min_row = (w(0)-1)/2;
  int max_row = nr - ((w(0)-1)/2);
  int min_col = (w(1)-1)/2;
  int max_col = nc - ((w(1)-1)/2);
  int n_elem = nr * nc;
  //int center_idx = ((w[0]*w[1])-1)/2;
  
  //Z = aX2 + bY2 + cXY + dX + eY + f
  
  NumericMatrix out = NumericMatrix(n_elem, 7);
  colnames(out)= CharacterVector::create("a", "b", "c", "d", "e", "f", "mask");
  out.fill(NA_REAL);
  out(_,6)=rep(0,100); //initialize mask with 0's
  
  //NEED AT LEAST 6 POINTS TO CALCULATE BECAUSE NEED AS MANY POINTS AS PARAMETERS
  int thresh = 6;
  for(int i = min_row; i< max_row; ++i) {
    for(int j = min_col; j < max_col; ++j){
      IntegerVector idx = IntegerVector(2);
      idx(0)= i;
      idx(1)=j;
      NumericMatrix curr_window = C_extract_window(r, w, idx);
      NumericVector Z = as<NumericVector>(curr_window);
      LogicalVector NA_idx = is_na(Z);
      int n_obs = sum(!NA_idx);
      
      if((is_true(any(NA_idx)) && (!na_rm)) || (n_obs < thresh)) {} else {
        NumericVector Z_trim_vect = Z[!NA_idx];
        NumericMatrix Z_trim(n_obs,1, Z_trim_vect.begin());
        NumericMatrix X_trim = subset_mat_rows(X, !NA_idx);
        
        int curr_elem_idx = i*nc + j; //rasters are indexed moving across rows
        NumericVector uni_Zvals = unique(Z_trim_vect);
        if(uni_Zvals.length() == 1){
          //If all Z values are the same, intercept should just be the value and all other parameters are 0. mask is 1 indicating all values are the same
          out(curr_elem_idx, 5) = uni_Zvals[0]; //f
          out(curr_elem_idx, 0) = 0; //a
          out(curr_elem_idx, 1) = 0; //b
          out(curr_elem_idx, 2) = 0; //c
          out(curr_elem_idx, 3) = 0; //d
          out(curr_elem_idx, 4) = 0; //e
          out(curr_elem_idx,6) = 1; //mask
          } else{
            NumericVector params = C_OLS_params(as<arma::mat>(X_trim), as<arma::mat>(Z_trim));
            out(curr_elem_idx, 5) =  params[0]; //f
            out(curr_elem_idx, 0) =  params[1]; //a
            out(curr_elem_idx, 1) =  params[2]; //b
            out(curr_elem_idx, 2) =  params[3]; //c
            out(curr_elem_idx, 3) =  params[4]; //d
            out(curr_elem_idx, 4) =  params[5]; //e
            }
          }}
  }
  return(out);
}

//Multiscale metrics across matrix using sliding window (Planar Fit SD)
// [[Rcpp::export]]
NumericVector C_multiscale1(NumericMatrix r, IntegerVector w, NumericMatrix X, bool na_rm){
  int nr= r.nrow();
  int nc= r.ncol();
  int min_row = (w(0)-1)/2;
  int max_row = nr - ((w(0)-1)/2);
  int min_col = (w(1)-1)/2;
  int max_col = nc - ((w(1)-1)/2);
  int n_elem = nr * nc;
  //int center_idx = ((w[0]*w[1])-1)/2;
  
  //Z = cXY + dX + eY + f
  NumericVector out = NumericVector(n_elem); //SD Residuals
  out.fill(NA_REAL);
  
  //NEED AT LEAST 3 POINTS TO CALCULATE BECAUSE NEED AS MANY POINTS AS PARAMETERS
  //SET THRESH TO 4 B/C WITH 3 RESIDUALS WILL ALWAYS BE 0
  
  int thresh = 4;
  for(int i = min_row; i< max_row; ++i) {
    for(int j = min_col; j < max_col; ++j){
      IntegerVector idx = IntegerVector(2);
      idx(0)= i;
      idx(1)=j;
      NumericMatrix curr_window = C_extract_window(r, w, idx);
      NumericVector Z = as<NumericVector>(curr_window);
      LogicalVector NA_idx = is_na(Z);
      int n_obs = sum(!NA_idx);
      
      if((is_true(any(NA_idx)) && (!na_rm)) || (n_obs < thresh)) {} else {
        NumericVector Z_trim_vect = Z[!NA_idx];
        NumericMatrix Z_trim(n_obs,1, Z_trim_vect.begin());
        NumericMatrix X_trim = subset_mat_rows(X, !NA_idx);
        
        int curr_elem_idx = i*nc + j; //rasters are indexed moving across rows
        NumericVector uni_Zvals = unique(Z_trim_vect);
        if(uni_Zvals.length() == 1){
          out[curr_elem_idx] = 0; //SD resid
          } else{
            NumericVector resid = C_OLS_resid(as<arma::mat>(X_trim), as<arma::mat>(Z_trim));
            out[curr_elem_idx] =  sd(resid); //SD resid
            }
          }}
  }
  return(out);
}


//FOR TERRA
// [[Rcpp::export]]
NumericMatrix C_Multiscale2b(NumericVector z, NumericMatrix X_full, bool na_rm, size_t ni, size_t nw) {
  
  size_t nlyr = X_full.ncol() + 1; //Number of layers
  size_t n = nlyr * ni;	//Number of elements in all layers
  NumericMatrix out = NumericMatrix(ni, 7);
  out.fill(NA_REAL);
  out(_,6)=rep(0,out.nrow()); //initialize mask with 0's
  colnames(out)= CharacterVector::create("a", "b", "c", "d", "e", "f", "mask");
  
  int thresh = 6; //NEED AT LEAST 6 POINTS TO CALCULATE BECAUSE NEED AS MANY POINTS AS PARAMETERS
  
  for (size_t i=0; i<ni; i++) {
    size_t start = i*nw;
    size_t end = start+nw-1;
    NumericVector zw_full = z[Rcpp::Range(start,end)]; //Current window of elevation values
    LogicalVector NA_idx = is_na(zw_full);
    int n_obs = sum(!NA_idx);
    if((is_true(any(NA_idx)) && (!na_rm)) || (n_obs < thresh)) {} else {
      NumericVector zw = zw_full[!NA_idx];
      NumericMatrix Z(n_obs,1, zw.begin());
      NumericMatrix X = subset_mat_rows(X_full, !NA_idx);
      NumericVector uni_Zvals = unique(zw);
      if(uni_Zvals.length() == 1){
        //If all Z values are the same, intercept should just be the value and all other parameters are 0. mask is 1 indicating all values are the same
        out(i, 5) = uni_Zvals[0]; //f
        out(i, 0) = 0; //a
        out(i, 1) = 0; //b
        out(i, 2) = 0; //c
        out(i, 3) = 0; //d
        out(i, 4) = 0; //e
        out(i,6) = 1; //mask
      } else{
        NumericVector params = C_OLS_params(as<arma::mat>(X), as<arma::mat>(Z));
        out(i, 5) =  params[0]; //f
        out(i, 0) =  params[1]; //a
        out(i, 1) =  params[2]; //b
        out(i, 2) =  params[3]; //c
        out(i, 3) =  params[4]; //d
        out(i, 4) =  params[5]; //e
      }
    }}
  return out;
}


// [[Rcpp::export]]
NumericMatrix C_OLSraw2(NumericVector y, arma::mat X, size_t ni, size_t nw) {
  
  arma::mat Xt = trans(X);
  arma::mat XtX = Xt * X;
  double d = det(XtX);
  
  size_t nlyr = X.n_cols;
  size_t n = nlyr * ni;
  NumericMatrix out(ni, nlyr);
  out.fill(NA_REAL);
  
  if(d==0){
    return out;
  } 
  arma::mat XtX_inv= inv(XtX);
  arma::mat H = X * XtX_inv * Xt;
  arma::mat Y(nw, 1);
  
  for (size_t i=0; i<ni; i++) {
    size_t start = i*nw;
    size_t end = start+nw-1;
    NumericVector yw = y[Rcpp::Range(start,end)];
    if (all(! Rcpp::is_na(yw))) {
      for(size_t j=0; j<nw; j++){
        Y(j,0) = yw[j];
      }
      NumericVector B = wrap(XtX_inv * (Xt * Y));	
      // output values must be interleaved by band (layer), not by pixel
      for (size_t j=0; j<nlyr; j++) { 
        out(i,j) = B[j];
      }
    }
  }
  
  return out;
}