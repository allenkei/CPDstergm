// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CPD_STERGM_cpp
List CPD_STERGM_cpp(int ADMM_iter, int theta_iter, int z_iter, List H_pos_list, List H_neg_list, List y_pos_list, List y_neg_list, arma::mat theta_mat, arma::mat z_mat, arma::mat u_mat, arma::mat X_mat, arma::vec d_vec, double alpha, int tau, int p1, int p2, double lambda, double theta_tol, double ADMM_tol, bool verbose);
RcppExport SEXP _CPDstergm_CPD_STERGM_cpp(SEXP ADMM_iterSEXP, SEXP theta_iterSEXP, SEXP z_iterSEXP, SEXP H_pos_listSEXP, SEXP H_neg_listSEXP, SEXP y_pos_listSEXP, SEXP y_neg_listSEXP, SEXP theta_matSEXP, SEXP z_matSEXP, SEXP u_matSEXP, SEXP X_matSEXP, SEXP d_vecSEXP, SEXP alphaSEXP, SEXP tauSEXP, SEXP p1SEXP, SEXP p2SEXP, SEXP lambdaSEXP, SEXP theta_tolSEXP, SEXP ADMM_tolSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type ADMM_iter(ADMM_iterSEXP);
    Rcpp::traits::input_parameter< int >::type theta_iter(theta_iterSEXP);
    Rcpp::traits::input_parameter< int >::type z_iter(z_iterSEXP);
    Rcpp::traits::input_parameter< List >::type H_pos_list(H_pos_listSEXP);
    Rcpp::traits::input_parameter< List >::type H_neg_list(H_neg_listSEXP);
    Rcpp::traits::input_parameter< List >::type y_pos_list(y_pos_listSEXP);
    Rcpp::traits::input_parameter< List >::type y_neg_list(y_neg_listSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta_mat(theta_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z_mat(z_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type u_mat(u_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_mat(X_matSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type d_vec(d_vecSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< int >::type p2(p2SEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type theta_tol(theta_tolSEXP);
    Rcpp::traits::input_parameter< double >::type ADMM_tol(ADMM_tolSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(CPD_STERGM_cpp(ADMM_iter, theta_iter, z_iter, H_pos_list, H_neg_list, y_pos_list, y_neg_list, theta_mat, z_mat, u_mat, X_mat, d_vec, alpha, tau, p1, p2, lambda, theta_tol, ADMM_tol, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CPDstergm_CPD_STERGM_cpp", (DL_FUNC) &_CPDstergm_CPD_STERGM_cpp, 20},
    {NULL, NULL, 0}
};

RcppExport void R_init_CPDstergm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}