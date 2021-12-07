#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


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
    NumericVector B = Rcpp::as<Rcpp::NumericVector>(wrap(XtX_inv * (Xt * Y)));
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

// [[Rcpp::export]]
NumericMatrix C_WoodEvans(NumericVector z, NumericMatrix X_full, bool na_rm, size_t ni, size_t nw) {
  
  size_t nlyr = X_full.ncol() + 1; //Number of layers
  NumericMatrix out = NumericMatrix(ni, nlyr);
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
        out(i, 0) = 0; //a
        out(i, 1) = 0; //b
        out(i, 2) = 0; //c
        out(i, 3) = 0; //d
        out(i, 4) = 0; //e
        out(i, 5) = uni_Zvals[0]; //f
        out(i, 6) = 1; //mask
      } else{
        NumericVector params = C_OLS_params(as<arma::mat>(X), as<arma::mat>(Z));
        out(i, 0) =  params[0]; //a
        out(i, 1) =  params[1]; //b
        out(i, 2) =  params[2]; //c
        out(i, 3) =  params[3]; //d
        out(i, 4) =  params[4]; //e
        out(i, 5) =  params[5]; //f
      }
    }}
  return out;
}

//Multiscale metrics across matrix using sliding window (Planar Fit SD)
// [[Rcpp::export]]
NumericVector C_AdjSD(NumericVector z, NumericMatrix X_full, bool na_rm, size_t ni, size_t nw){
  NumericVector out = NumericVector(ni, NA_REAL);
  //Z = dX + eY + f
  
  //NEED AT LEAST 3 POINTS TO CALCULATE BECAUSE NEED AS MANY POINTS AS PARAMETERS
  //SET THRESH TO 4 B/C WITH 3 RESIDUALS WILL ALWAYS BE 0
  
  int thresh = 4;
  for (size_t i=0; i< ni; i++) {
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
        out[i] = 0; //If all Z values are the same residuals are 0
        } else{
          NumericVector resid = C_OLS_resid(as<arma::mat>(X), as<arma::mat>(Z));
          out[i] =  sd(resid); //SD resid
          }
    }}
  return out;
}

// Calculate area of triangle based on side lengths
// [[Rcpp::export]]
double C_TriArea (double a, double b, double c){
  double s = (a+b+c)/2;
  double out =sqrt(s*(s-a)*(s-b)*(s-c));
  return out;
}

//Surface Area
// [[Rcpp::export]]
NumericVector C_SurfaceArea (NumericVector z, double x_res, double y_res, size_t ni, size_t nw){
  NumericVector out = NumericVector(ni, NA_REAL);
  double Lx2= pow(x_res, 2);
  double Ly2= pow(y_res, 2);
  double Ld2= Lx2 + Ly2;
  
  for (size_t i=0; i< ni; i++) {
    size_t start = i*nw;
    size_t end = start+nw-1;
    NumericVector zw = z[Rcpp::Range(start,end)]; //Current window of elevation values
    //|A|B|C|
    //|D|E|F|
    //|G|H|I|
    //Calculate Edge Lengths
    //Horiz
    double AB= sqrt(Lx2+pow(zw[0]-zw[1],2))/2;
    double BC= sqrt(Lx2+pow(zw[1]-zw[2],2))/2;
    double DE= sqrt(Lx2+pow(zw[3]-zw[4],2))/2;
    double EF= sqrt(Lx2+pow(zw[4]-zw[5],2))/2;
    double GH= sqrt(Lx2+pow(zw[6]-zw[7],2))/2;
    double HI= sqrt(Lx2+pow(zw[7]-zw[8],2))/2;
    //Vertical
    double AD= sqrt(Ly2+pow(zw[0]-zw[3],2))/2;
    double BE= sqrt(Ly2+pow(zw[1]-zw[4],2))/2;
    double CF= sqrt(Ly2+pow(zw[2]-zw[5],2))/2;
    double DG= sqrt(Ly2+pow(zw[3]-zw[6],2))/2;
    double EH= sqrt(Ly2+pow(zw[4]-zw[7],2))/2;
    double FI= sqrt(Ly2+pow(zw[5]-zw[8],2))/2;
    //Diagonal
    double EA= sqrt(Ld2+pow(zw[4]-zw[0],2))/2;
    double EC= sqrt(Ld2+pow(zw[4]-zw[2],2))/2;
    double EG= sqrt(Ld2+pow(zw[4]-zw[6],2))/2;
    double EI= sqrt(Ld2+pow(zw[4]-zw[8],2))/2;
    
    //Sum area of the 8 triangles
    out[i] = C_TriArea(EA, AB, BE) + C_TriArea(BE, BC, EC) + C_TriArea(AD, DE, EA) + C_TriArea(EC, CF, EF) + C_TriArea(DE, DG, EG) + C_TriArea(EF, FI, EI) + C_TriArea(EG, EH, GH) + C_TriArea(EH, EI, HI);
  }
  return out;
}

//Count values
// [[Rcpp::export]]
NumericVector C_CountVals (NumericVector z, size_t ni, size_t nw){
  NumericVector out = NumericVector(ni, NA_REAL);
  for (size_t i=0; i< ni; i++) {
    size_t start = i*nw;
    size_t end = start+nw-1;
    NumericVector zw = z[Rcpp::Range(start,end)]; //Current window of elevation values
    out[i] = sum(!is_na(zw));
  }
  return out;
}