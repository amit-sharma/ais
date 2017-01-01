#include <Rcpp.h>
using namespace Rcpp;

typedef NumericVector (*rDistrFnPtr)(NumericVector x , int j, double mult);
typedef double (*dCondDensityFnPtr)(NumericVector x, NumericVector y, int n, double mult);