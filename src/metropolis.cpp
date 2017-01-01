#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "probability_functions.h"
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
//RNGScope scope;
//IntegerVector powers_dirichlet= IntegerVector::create(1, 2, 4, 5, 6, 1, 5, 8);

// [[Rcpp::export]]
NumericVector metropolisC(NumericVector x, double beta, int num_iterations_mcmc, 
                          NumericVector other_params) {
  int n = x.size();
  NumericVector proposal(n);
  double log_old_p, log_new_p;
  double new_fa, log_diff;
  NumericVector rand_num;
  for(int i =1; i<=num_iterations_mcmc; i++){
    for(int j=1; j<=9; j=j*3){
        proposal = x + rnorm(n, 0, 0.05*j);
        new_fa = faC(proposal);
        if(new_fa != -INFINITY){
          //std::cout<<x<<"ggf"<<proposal<<std::endl;
          //std::cout<<faC(x)<<new_fa<<std::endl;
          //std::cout<<pow(exp(faC(x)),(1-beta));
          log_old_p = faC(x)*(1-beta) + fbC(x, other_params)*beta;
          log_new_p = new_fa*(1-beta) + fbC(proposal, other_params)*beta;
          log_diff = log_new_p - log_old_p;
          //std::cout<<log_new_p-log_old_p<<" "<<log_new_p<<log_old_p<<std::endl;
          if(log_new_p != -INFINITY){
            if(exp(log_diff) > runif(1)[0]){
              for(int index=0; index<n; index++){
                x[index] = proposal[index];
              }
            }
          }
        }
        
    }
  }
  return x;
}

// [[Rcpp::export]]
NumericVector metropolisC2(NumericVector x, double beta, int num_iterations_mcmc,
                           NumericVector other_params) {
  int n = x.size();
  NumericVector proposal(n);
  double log_old_p, log_new_p, log_diff;
  double new_fa;
  bool do_compute_old_p=true;
  int num_potential_changes=0;
  for(int i =1; i<=num_iterations_mcmc; i++){
    for(int j=1; j<=9; j=j*3){
      proposal = x + rnorm(n, 0, 0.05*j);
      //std::cout<<x<<"ggf"<<proposal<<std::endl;
      //std::cout<<pow(exp(faC(x)),(1-beta));
      new_fa = faC(proposal);
      if(new_fa != -INFINITY){
        if(do_compute_old_p){
          log_old_p = faC(x)*(1-beta) + fbC(x, other_params)*beta;
          do_compute_old_p = false;
        }
        log_new_p = new_fa*(1-beta) + fbC(proposal, other_params)*beta;
        
        log_diff = log_new_p - log_old_p;
        //std::cout<<new_p/old_p<<std::endl; 
        if(log_new_p != -INFINITY){
          num_potential_changes++;
          if(exp(log_diff) > runif(1)[0]){
            for(int index=0; index<n; index++){
              x[index] = proposal[index];
            }
            do_compute_old_p=true;
            
          }
        }
      }
      
    }
  }
  //std::cout<<num_potential_changes<<" Number of potential changes"<<std::endl;
  return x;
}


double dbeta_vectorized(NumericVector x, NumericVector a, NumericVector b, bool logscale){
  //logscale is non-zero for outputting log
  double log_g = 0;
  for(int index=0;index<x.size();index++){
    log_g = log_g + R::dbeta(x[index], a[index], b[index], logscale);
  }
  return log_g;
  
}
// [[Rcpp::export]]
NumericVector metropolisCbeta(NumericVector x, double beta, int num_iterations_mcmc,
                           NumericVector other_params) {
  int n = x.size();
  NumericVector proposal(n);
  double log_old_p, log_new_p, log_diff;
  double log_old_g, log_new_g;
  double new_fa;
  bool do_compute_old_p=true;
  int num_potential_changes=0;
  for(int i =1; i<=num_iterations_mcmc; i++){
    for(int j=1000; j<=100000; j=j*10){
      for(int index=0;index<n;index++){
        proposal[index] = rbeta(1, j*x[index], j*(1-x[index]))[0];
      }
      //std::cout<<x<<"ggf"<<proposal<<std::endl;
      //std::cout<<pow(exp(faC(x)),(1-beta));
      new_fa = faC(proposal);
      if(new_fa != -INFINITY){
        if(do_compute_old_p){
          log_old_p = faC(x)*(1-beta) + fbC(x, other_params)*beta;
          do_compute_old_p = false;
        }
        log_new_p = new_fa*(1-beta) + fbC(proposal, other_params)*beta;
        log_new_g = dbeta_vectorized(proposal, j*x, j*(1-x), true);
        log_old_g = dbeta_vectorized(x, j*proposal, j*(1-proposal), true);
        log_diff = (log_new_p - log_old_p) + (log_old_g -log_new_g);
        //std::cout<<new_p/old_p<<std::endl; 
        if(log_new_p != -INFINITY){
          num_potential_changes++;
          if(exp(log_diff) > runif(1)[0]){
            for(int index=0; index<n; index++){
              x[index] = proposal[index];
            }
            do_compute_old_p=true;
            
          }
        } /*else {
          std::cout<<log_new_p<<new_fa<<" "<<fbC(proposal, other_params)<<log_old_p<<"  gf  "<<proposal<<std::endl;
        }*/
      }
      
    }
  }
  //std::cout<<num_potential_changes<<" Number of potential changes"<<std::endl;
  return x;
}


/*** R code
log_likelihood <- function(theta) {
  log_like = 0
  if(DISTR=="dirichlet"){
    log_like = powers_dirichlet %*% log(theta)
    normalizing_constant = gamma(sum(powers_dirichlet+1))/prod(gamma(powers_dirichlet+1))
  } else if(DISTR=="hyperdirichlet") {
    theta_mat = matrix(theta, ncol=4,byrow=TRUE)
    theta_sum_terms = rowSums(theta_mat)
    log_like=powers_dirichlet %*% log(theta_sum_terms)
    normalizing_constant = (gamma(sum(powers_dirichlet+4)) * gamma(4)^8)/prod(gamma(powers_dirichlet+4))
  }
  
  log_like = log_like + log(normalizing_constant)
    return(log_like)
} 

fb <- function(e){
  e <- c(1,e)
  p <- e_to_p(e)
  if(any(p<0) || any(p>1) || sum(p)>1)
    return(log(0))
    
    out <- log_likelihood(p)
    return(log(Jacobian(e))+ out)
}




*/

/*** R code
powers_dirichlet= c(1, 2, 4, 5, 6, 1, 5, 8)
if(DISTR=='dirichlet'){
 n=length(powers_dirichlet)-1 #            # Number of dimensions of distribution (1 less than 10)
} else if (DISTR=="hyperdirichlet"){
 n=length(powers_dirichlet)*4 -1
} else {
 print("Error in distribution")
}
 
e = runif(n)

Jacobian(c(1,e))
JacobianC(e)
#log_likelihood(e_to_p(c(1,e)))
fa(e)
faC(e)
fb(e)
fbC(e)
# hack in this line: check length of rnorm
set.seed(1)
metropolis(e, fa=fa, fb=fb, beta=0.2, jump=function(x){rnorm(length(powers_dirichlet)-1, mean=0, sd=0.05*x)}, num_iterations_mcmc=10000)
set.seed(1)
metropolisC(e, 0.2, 10000) 
*/