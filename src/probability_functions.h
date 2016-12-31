#ifndef PROB_FUNC_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define PROB_FUNC_H

#include <Rcpp.h>
using namespace Rcpp;

double faC(NumericVector e);
double fbC(NumericVector e, NumericVector other_params);

#endif