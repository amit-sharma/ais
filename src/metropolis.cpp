#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

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
double faC(NumericVector e){
 if(is_true(any(e<0)))
   return -INFINITY;
 if(is_true(any(e>1)))
   return -INFINITY;
  return 0;
   
}
// [[Rcpp::export]]
double JacobianC(NumericVector e){
  double res=1;
  int len_vec = e.size();
  for(int i=0; i<len_vec-1; i++){
    res = res * pow(1-e[i], len_vec-i-1);
  }
  return(res);
}

// [[Rcpp::export]]
double log_likelihoodC(NumericVector theta, IntegerVector powers_dirichlet){
  double log_like, log_normalizing_constant;
  log_like = sum(as<NumericVector>(powers_dirichlet) * log(theta));
  log_normalizing_constant = lgamma(sum(powers_dirichlet+1)) - sum(lgamma(powers_dirichlet+1));
   
  // For hyperdirichlet
  /*NumericVector theta_sum_terms(theta.size()/4);
  for(int index=0;index < theta.size();index= index + 4){
    theta_sum_terms[index/4] = theta[index] + theta[index+1] + theta[index+2] + theta[index+3];
  }
  
  log_like=sum(as<NumericVector>(powers_dirichlet) * log(theta_sum_terms));
  //NumericVector temp=log(theta_sum_terms);std::cout<<temp<<std::endl;
  log_normalizing_constant = lgamma(sum(powers_dirichlet+4)) + 8*lgamma(4) - sum(lgamma(powers_dirichlet+4));
   */
  if(!(log_like >-90)){
    //std::cout<<log_like<<powers_dirichlet<<std::endl;
    //std::cout<<theta_sum_terms<<theta<<std::endl;
  }
  log_like = log_like + log_normalizing_constant;
  return log_like;
}

// [[Rcpp::export]]
NumericVector e_to_pC(NumericVector e){
  int len_vec = e.size()+1;
  NumericVector p(len_vec);
  double mult=1;
  
  for(int i=0; i<len_vec-1; i++){
    p[i] = e[i]*mult;
    mult = mult*(1-e[i]);
  }
  p[len_vec-1]=mult;
  return p;
}

// [[Rcpp::export]]
double fbC(NumericVector e, NumericVector other_params){
  //Already checked in faC
 if(is_true(any(e<0)))
    return -INFINITY ;
 if(is_true(any(e>1)))
      return -INFINITY;
  double jacobian, out;
  IntegerVector powers_dirichlet= as<IntegerVector>(other_params); // because we know that powers will be integers
  NumericVector p;
  p = e_to_pC(e);
  
  if(is_true(any(p<0)))
    std::cout<<"ttr"<<e<<p<<std::endl;
  out = log_likelihoodC(p, powers_dirichlet);
  jacobian = JacobianC(e);
  //std::cout<<powers_dirichlet<<std::endl;
  //std::cout<<log(jacobian)<<" "<<out<<std::endl;
  return log(jacobian) + out;
}

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