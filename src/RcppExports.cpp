// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "ais_types.h"
#include "../inst/include/ais.h"
#include <Rcpp.h>

using namespace Rcpp;

// metropolisC
NumericVector metropolisC(NumericVector x, double beta, int num_iterations_mcmc, List other_params);
RcppExport SEXP ais_metropolisC(SEXP xSEXP, SEXP betaSEXP, SEXP num_iterations_mcmcSEXP, SEXP other_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type num_iterations_mcmc(num_iterations_mcmcSEXP);
    Rcpp::traits::input_parameter< List >::type other_params(other_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolisC(x, beta, num_iterations_mcmc, other_params));
    return rcpp_result_gen;
END_RCPP
}
// metropolisC2
List metropolisC2(NumericVector x, double beta, int num_iterations_mcmc, SEXP rproposal_fn_xpsexp, SEXP dproposal_fn_xpsexp, List other_params);
RcppExport SEXP ais_metropolisC2(SEXP xSEXP, SEXP betaSEXP, SEXP num_iterations_mcmcSEXP, SEXP rproposal_fn_xpsexpSEXP, SEXP dproposal_fn_xpsexpSEXP, SEXP other_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type num_iterations_mcmc(num_iterations_mcmcSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rproposal_fn_xpsexp(rproposal_fn_xpsexpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dproposal_fn_xpsexp(dproposal_fn_xpsexpSEXP);
    Rcpp::traits::input_parameter< List >::type other_params(other_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolisC2(x, beta, num_iterations_mcmc, rproposal_fn_xpsexp, dproposal_fn_xpsexp, other_params));
    return rcpp_result_gen;
END_RCPP
}
// metropolisCbeta
List metropolisCbeta(NumericVector x, double beta, int num_iterations_mcmc, SEXP rproposal_fn_xpsexp, SEXP dproposal_fn_xpsexp, List other_params, int debug);
RcppExport SEXP ais_metropolisCbeta(SEXP xSEXP, SEXP betaSEXP, SEXP num_iterations_mcmcSEXP, SEXP rproposal_fn_xpsexpSEXP, SEXP dproposal_fn_xpsexpSEXP, SEXP other_paramsSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type num_iterations_mcmc(num_iterations_mcmcSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rproposal_fn_xpsexp(rproposal_fn_xpsexpSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dproposal_fn_xpsexp(dproposal_fn_xpsexpSEXP);
    Rcpp::traits::input_parameter< List >::type other_params(other_paramsSEXP);
    Rcpp::traits::input_parameter< int >::type debug(debugSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolisCbeta(x, beta, num_iterations_mcmc, rproposal_fn_xpsexp, dproposal_fn_xpsexp, other_params, debug));
    return rcpp_result_gen;
END_RCPP
}
// faC
double faC(NumericVector e);
RcppExport SEXP ais_faC(SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(faC(e));
    return rcpp_result_gen;
END_RCPP
}
// JacobianC
double JacobianC(NumericVector e);
RcppExport SEXP ais_JacobianC(SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(JacobianC(e));
    return rcpp_result_gen;
END_RCPP
}
// logJacobianC
double logJacobianC(NumericVector e);
RcppExport SEXP ais_logJacobianC(SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(logJacobianC(e));
    return rcpp_result_gen;
END_RCPP
}
// log_likelihoodC
double log_likelihoodC(NumericVector theta, IntegerVector powers_dirichlet, Rcpp::List theta_sum_list);
RcppExport SEXP ais_log_likelihoodC(SEXP thetaSEXP, SEXP powers_dirichletSEXP, SEXP theta_sum_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type powers_dirichlet(powers_dirichletSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type theta_sum_list(theta_sum_listSEXP);
    rcpp_result_gen = Rcpp::wrap(log_likelihoodC(theta, powers_dirichlet, theta_sum_list));
    return rcpp_result_gen;
END_RCPP
}
// e_to_pC
NumericVector e_to_pC(NumericVector e);
RcppExport SEXP ais_e_to_pC(SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(e_to_pC(e));
    return rcpp_result_gen;
END_RCPP
}
// fbC
double fbC(NumericVector e, List other_params);
RcppExport SEXP ais_fbC(SEXP eSEXP, SEXP other_paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type e(eSEXP);
    Rcpp::traits::input_parameter< List >::type other_params(other_paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(fbC(e, other_params));
    return rcpp_result_gen;
END_RCPP
}
// rproposal_distr
NumericVector rproposal_distr(NumericVector x, int n, double mult);
RcppExport SEXP ais_rproposal_distr(SEXP xSEXP, SEXP nSEXP, SEXP multSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mult(multSEXP);
    rcpp_result_gen = Rcpp::wrap(rproposal_distr(x, n, mult));
    return rcpp_result_gen;
END_RCPP
}
// dbeta_cond_distr
double dbeta_cond_distr(NumericVector x, NumericVector y, int n, double mult);
RcppExport SEXP ais_dbeta_cond_distr(SEXP xSEXP, SEXP ySEXP, SEXP nSEXP, SEXP multSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mult(multSEXP);
    rcpp_result_gen = Rcpp::wrap(dbeta_cond_distr(x, y, n, mult));
    return rcpp_result_gen;
END_RCPP
}
// putFunPtrInXPtr
XPtr<rDistrFnPtr> putFunPtrInXPtr(std::string distr_name);
RcppExport SEXP ais_putFunPtrInXPtr(SEXP distr_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type distr_name(distr_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(putFunPtrInXPtr(distr_name));
    return rcpp_result_gen;
END_RCPP
}
// getCondDensityFuncXPtr
XPtr<dCondDensityFnPtr> getCondDensityFuncXPtr(std::string distr_name);
RcppExport SEXP ais_getCondDensityFuncXPtr(SEXP distr_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type distr_name(distr_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(getCondDensityFuncXPtr(distr_name));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP ais_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"ais_metropolisC", (DL_FUNC) &ais_metropolisC, 4},
    {"ais_metropolisC2", (DL_FUNC) &ais_metropolisC2, 6},
    {"ais_metropolisCbeta", (DL_FUNC) &ais_metropolisCbeta, 7},
    {"ais_faC", (DL_FUNC) &ais_faC, 1},
    {"ais_JacobianC", (DL_FUNC) &ais_JacobianC, 1},
    {"ais_logJacobianC", (DL_FUNC) &ais_logJacobianC, 1},
    {"ais_log_likelihoodC", (DL_FUNC) &ais_log_likelihoodC, 3},
    {"ais_e_to_pC", (DL_FUNC) &ais_e_to_pC, 1},
    {"ais_fbC", (DL_FUNC) &ais_fbC, 2},
    {"ais_rproposal_distr", (DL_FUNC) &ais_rproposal_distr, 3},
    {"ais_dbeta_cond_distr", (DL_FUNC) &ais_dbeta_cond_distr, 4},
    {"ais_putFunPtrInXPtr", (DL_FUNC) &ais_putFunPtrInXPtr, 1},
    {"ais_getCondDensityFuncXPtr", (DL_FUNC) &ais_getCondDensityFuncXPtr, 1},
    {"ais_rcpp_hello_world", (DL_FUNC) &ais_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_ais(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
