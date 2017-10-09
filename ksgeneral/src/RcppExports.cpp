// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cont_ks_distribution_Rcpp
double cont_ks_distribution_Rcpp(double n);
RcppExport SEXP _ksgeneral_cont_ks_distribution_Rcpp(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(cont_ks_distribution_Rcpp(n));
    return rcpp_result_gen;
END_RCPP
}
// cont_ks_distribution_Rcpp_alternative
double cont_ks_distribution_Rcpp_alternative(double n, double x);
RcppExport SEXP _ksgeneral_cont_ks_distribution_Rcpp_alternative(SEXP nSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(cont_ks_distribution_Rcpp_alternative(n, x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ksgeneral_cont_ks_distribution_Rcpp", (DL_FUNC) &_ksgeneral_cont_ks_distribution_Rcpp, 1},
    {"_ksgeneral_cont_ks_distribution_Rcpp_alternative", (DL_FUNC) &_ksgeneral_cont_ks_distribution_Rcpp_alternative, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ksgeneral(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}