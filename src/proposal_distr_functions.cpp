#include <Rcpp.h>
using namespace Rcpp;

#include "ais_types.h"



//[[Rcpp::export]]
NumericVector rproposal_distr(NumericVector x, int n, double mult){
  NumericVector params1 = mult*x;
  NumericVector params2 = mult*(1-x);
  NumericVector proposal(n);
  for(int index=0;index<n;index++){
    proposal[index] = R::rbeta(params1[index],params2[index]);
  }
  return proposal;
  
}

//[[Rcpp::export]]
double dbeta_cond_distr(NumericVector x, NumericVector y, int n, double mult){
  //logscale is non-zero for outputting log
  double log_g = 0;
  NumericVector a = mult*y;
  NumericVector b = mult*(1-y);
  
  for(int index=0;index<n;index++){
    log_g = log_g + R::dbeta(x[index], a[index], b[index], 1);
  }
  return log_g;
  
}


//[[Rcpp::export]]
XPtr<rDistrFnPtr> putFunPtrInXPtr(std::string distr_name) {
  if(distr_name=="default"){
    return(XPtr<rDistrFnPtr>(new rDistrFnPtr(&rproposal_distr)));
  } else {
    return XPtr<rDistrFnPtr>(R_NilValue); // runtime error as NULL no XPtr
  }
}


//[[Rcpp::export]]
XPtr<dCondDensityFnPtr> getCondDensityFuncXPtr(std::string distr_name) {
  if(distr_name=="default"){
    return(XPtr<dCondDensityFnPtr>(new dCondDensityFnPtr(&dbeta_cond_distr)));
  } else {
    return XPtr<dCondDensityFnPtr>(R_NilValue); // runtime error as NULL no XPtr
  }
}

