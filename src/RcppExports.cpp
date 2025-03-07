// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// subset_mat_rows
Rcpp::NumericMatrix subset_mat_rows(Rcpp::NumericMatrix Input_Matrix, Rcpp::LogicalVector Input_Log_Vec);
RcppExport SEXP _MultiscaleDTM_subset_mat_rows(SEXP Input_MatrixSEXP, SEXP Input_Log_VecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Input_Matrix(Input_MatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type Input_Log_Vec(Input_Log_VecSEXP);
    rcpp_result_gen = Rcpp::wrap(subset_mat_rows(Input_Matrix, Input_Log_Vec));
    return rcpp_result_gen;
END_RCPP
}
// C_OLS_params
NumericVector C_OLS_params(arma::mat X, arma::mat Y);
RcppExport SEXP _MultiscaleDTM_C_OLS_params(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(C_OLS_params(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// C_OLS_params2
NumericVector C_OLS_params2(arma::mat Xt, arma::mat XtX_inv, arma::mat Y);
RcppExport SEXP _MultiscaleDTM_C_OLS_params2(SEXP XtSEXP, SEXP XtX_invSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtX_inv(XtX_invSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(C_OLS_params2(Xt, XtX_inv, Y));
    return rcpp_result_gen;
END_RCPP
}
// C_OLS_resid
NumericVector C_OLS_resid(arma::mat X, arma::mat Y);
RcppExport SEXP _MultiscaleDTM_C_OLS_resid(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(C_OLS_resid(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// C_OLS_resid2
NumericVector C_OLS_resid2(arma::mat X, arma::mat Xt, arma::mat XtX_inv, arma::mat Y);
RcppExport SEXP _MultiscaleDTM_C_OLS_resid2(SEXP XSEXP, SEXP XtSEXP, SEXP XtX_invSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtX_inv(XtX_invSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(C_OLS_resid2(X, Xt, XtX_inv, Y));
    return rcpp_result_gen;
END_RCPP
}
// C_Qfit1_narmT
NumericMatrix C_Qfit1_narmT(NumericVector z, NumericMatrix X_full, size_t ni, size_t nw);
RcppExport SEXP _MultiscaleDTM_C_Qfit1_narmT(SEXP zSEXP, SEXP X_fullSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X_full(X_fullSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Qfit1_narmT(z, X_full, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// C_Qfit1_narmF
NumericMatrix C_Qfit1_narmF(NumericVector z, arma::mat X, arma::mat Xt, arma::mat XtX_inv, size_t ni, size_t nw);
RcppExport SEXP _MultiscaleDTM_C_Qfit1_narmF(SEXP zSEXP, SEXP XSEXP, SEXP XtSEXP, SEXP XtX_invSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtX_inv(XtX_invSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Qfit1_narmF(z, X, Xt, XtX_inv, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// C_Qfit2_narmT
NumericMatrix C_Qfit2_narmT(NumericVector z, NumericMatrix X_full, size_t ni, size_t nw);
RcppExport SEXP _MultiscaleDTM_C_Qfit2_narmT(SEXP zSEXP, SEXP X_fullSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X_full(X_fullSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Qfit2_narmT(z, X_full, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// C_Qfit2_narmF
NumericMatrix C_Qfit2_narmF(NumericVector z, arma::mat X, arma::mat Xt, arma::mat XtX_inv, size_t ni, size_t nw);
RcppExport SEXP _MultiscaleDTM_C_Qfit2_narmF(SEXP zSEXP, SEXP XSEXP, SEXP XtSEXP, SEXP XtX_invSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtX_inv(XtX_invSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Qfit2_narmF(z, X, Xt, XtX_inv, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// C_AdjSD_narmT
NumericVector C_AdjSD_narmT(NumericVector z, NumericMatrix X_full, size_t ni, size_t nw);
RcppExport SEXP _MultiscaleDTM_C_AdjSD_narmT(SEXP zSEXP, SEXP X_fullSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type X_full(X_fullSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_AdjSD_narmT(z, X_full, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// C_AdjSD_narmF
NumericVector C_AdjSD_narmF(NumericVector z, arma::mat X, arma::mat Xt, arma::mat XtX_inv, size_t ni, size_t nw);
RcppExport SEXP _MultiscaleDTM_C_AdjSD_narmF(SEXP zSEXP, SEXP XSEXP, SEXP XtSEXP, SEXP XtX_invSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtX_inv(XtX_invSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_AdjSD_narmF(z, X, Xt, XtX_inv, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// C_Pfit1_narmF
NumericMatrix C_Pfit1_narmF(NumericVector z, arma::mat X, arma::mat Xt, arma::mat XtX_inv, LogicalVector idx, size_t ni, size_t nw);
RcppExport SEXP _MultiscaleDTM_C_Pfit1_narmF(SEXP zSEXP, SEXP XSEXP, SEXP XtSEXP, SEXP XtX_invSEXP, SEXP idxSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xt(XtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type XtX_inv(XtX_invSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_Pfit1_narmF(z, X, Xt, XtX_inv, idx, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// C_TriArea
double C_TriArea(double a, double b, double c);
RcppExport SEXP _MultiscaleDTM_C_TriArea(SEXP aSEXP, SEXP bSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(C_TriArea(a, b, c));
    return rcpp_result_gen;
END_RCPP
}
// C_SurfaceArea
NumericVector C_SurfaceArea(NumericVector z, double x_res, double y_res, bool na_rm, size_t ni, size_t nw);
RcppExport SEXP _MultiscaleDTM_C_SurfaceArea(SEXP zSEXP, SEXP x_resSEXP, SEXP y_resSEXP, SEXP na_rmSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type x_res(x_resSEXP);
    Rcpp::traits::input_parameter< double >::type y_res(y_resSEXP);
    Rcpp::traits::input_parameter< bool >::type na_rm(na_rmSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_SurfaceArea(z, x_res, y_res, na_rm, ni, nw));
    return rcpp_result_gen;
END_RCPP
}
// C_CountVals
NumericVector C_CountVals(NumericVector z, size_t ni, size_t nw);
RcppExport SEXP _MultiscaleDTM_C_CountVals(SEXP zSEXP, SEXP niSEXP, SEXP nwSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< size_t >::type ni(niSEXP);
    Rcpp::traits::input_parameter< size_t >::type nw(nwSEXP);
    rcpp_result_gen = Rcpp::wrap(C_CountVals(z, ni, nw));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MultiscaleDTM_subset_mat_rows", (DL_FUNC) &_MultiscaleDTM_subset_mat_rows, 2},
    {"_MultiscaleDTM_C_OLS_params", (DL_FUNC) &_MultiscaleDTM_C_OLS_params, 2},
    {"_MultiscaleDTM_C_OLS_params2", (DL_FUNC) &_MultiscaleDTM_C_OLS_params2, 3},
    {"_MultiscaleDTM_C_OLS_resid", (DL_FUNC) &_MultiscaleDTM_C_OLS_resid, 2},
    {"_MultiscaleDTM_C_OLS_resid2", (DL_FUNC) &_MultiscaleDTM_C_OLS_resid2, 4},
    {"_MultiscaleDTM_C_Qfit1_narmT", (DL_FUNC) &_MultiscaleDTM_C_Qfit1_narmT, 4},
    {"_MultiscaleDTM_C_Qfit1_narmF", (DL_FUNC) &_MultiscaleDTM_C_Qfit1_narmF, 6},
    {"_MultiscaleDTM_C_Qfit2_narmT", (DL_FUNC) &_MultiscaleDTM_C_Qfit2_narmT, 4},
    {"_MultiscaleDTM_C_Qfit2_narmF", (DL_FUNC) &_MultiscaleDTM_C_Qfit2_narmF, 6},
    {"_MultiscaleDTM_C_AdjSD_narmT", (DL_FUNC) &_MultiscaleDTM_C_AdjSD_narmT, 4},
    {"_MultiscaleDTM_C_AdjSD_narmF", (DL_FUNC) &_MultiscaleDTM_C_AdjSD_narmF, 6},
    {"_MultiscaleDTM_C_Pfit1_narmF", (DL_FUNC) &_MultiscaleDTM_C_Pfit1_narmF, 7},
    {"_MultiscaleDTM_C_TriArea", (DL_FUNC) &_MultiscaleDTM_C_TriArea, 3},
    {"_MultiscaleDTM_C_SurfaceArea", (DL_FUNC) &_MultiscaleDTM_C_SurfaceArea, 6},
    {"_MultiscaleDTM_C_CountVals", (DL_FUNC) &_MultiscaleDTM_C_CountVals, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_MultiscaleDTM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
