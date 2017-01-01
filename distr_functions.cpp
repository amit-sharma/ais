//[[Rcpp::depends(ais)]]
#include <Rcpp.h>
using namespace Rcpp;

typedef NumericVector (*rDistrFnPtr)(NumericVector x , int j, double mult);
typedef double (*dCondDensityFnPtr)(NumericVector x, NumericVector y, int n, double mult);


//[[Rcpp::export]]
NumericVector rbeta_distr(NumericVector x, int n, double mult){
  NumericVector params1 = mult*x;
  NumericVector params2 = mult*(1-x);
  NumericVector proposal(n);
  for(int index=0;index<n;index++){
    proposal[index] = R::rbeta(params1[index],params2[index]);
  }
  return proposal;
  
}

//[[Rcpp::export]]
XPtr<rDistrFnPtr> putFunPtrInXPtrCustom() {
  return(XPtr<rDistrFnPtr>(new rDistrFnPtr(&rbeta_distr)));
}