#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

#include "ais_types.h"
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
//typedef NumericVector (*funcPtr)(NumericVector x , int j, double mult);

//Do not use this function. superseded by the C2 or Cbeta function below.
// [[Rcpp::export]]
NumericVector metropolisC(NumericVector x, double beta, int num_iterations_mcmc, 
                          List other_params) {
  int n = x.size();
  NumericVector proposal(n);
  double log_old_p, log_new_p;
  double new_fa, log_diff;
  NumericVector rand_num;
  int successes=0;
  NumericVector covariances = NumericVector::create(0.05,0.15,0.4, 0.7);
  for(int i =1; i<=num_iterations_mcmc; i++){
    //for(int j=1; j<=9; j=j*3){
    for(int j=0; j<covariances.size();j++){
        proposal = x + rnorm(n, 0, covariances[j]);//0.05*j);
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
            successes += 1;
            if(exp(log_diff) > runif(1)[0]){
              for(int index=0; index<n; index++){
                x[index] = proposal[index];
              }
            }
          }
        }
        
    }
  }
  //std::cout<<"Acceptance Ratio"<<successes/(num_iterations_mcmc*covariances.size());
  return x;
}

// [[Rcpp::export]]
List metropolisC2(NumericVector x, double beta, int num_iterations_mcmc,
                           SEXP rproposal_fn_xpsexp, SEXP dproposal_fn_xpsexp,
                           List other_params) {
  int n = x.size();
  NumericVector proposal(n);
  double log_old_p, log_new_p, log_diff;
  double new_fa;
  bool do_compute_old_p=true;
  int num_potential_changes=0;
  int num_actual_changes=0;
  int num_proposals_indomain = 0;
  int num_proposals = 0;
  NumericVector covariances = NumericVector::create(0.01,0.05,0.1, 0.3);
  for(int i =1; i<=num_iterations_mcmc; i++){
    for(int j=0; j<covariances.size(); j++){
      proposal = x + rnorm(n, 0, covariances[j]);//0.05*j);
      num_proposals++;
      //std::cout<<x<<"ggf"<<proposal<<std::endl;
      //std::cout<<pow(exp(faC(x)),(1-beta));
      new_fa = faC(proposal);
      if(new_fa != -INFINITY){
        num_proposals_indomain += 1;
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
            num_actual_changes++;
            for(int index=0; index<n; index++){
              x[index] = proposal[index];
            }
            do_compute_old_p=true;
            
          }
        }
      }
      
    }
  }
  //std::cout<<num_proposals_indomain<<" "<<num_potential_changes<<" Number of potential changes "<<1.0*num_potential_changes/(num_iterations_mcmc*covariances.size())<<" "<<1.0*num_proposals_indomain/(num_iterations_mcmc*covariances.size())<<" "<<1.0*num_actual_changes/(num_iterations_mcmc*covariances.size())<<" "<<1.0*num_actual_changes/num_potential_changes<<std::endl;
  return List::create(
    _["x_vec"]=x,
    _["num_proposals"]=num_proposals,
    _["num_proposals_indomain"]=num_proposals_indomain,
    _["num_potential_changes"]=num_potential_changes,
    _["num_actual_changes"]=num_actual_changes
  );
}


// From http://gallery.rcpp.org/articles/passing-cpp-function-pointers/

/*
// [[Rcpp::export]]
XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
  if (fstr == "fun1")
    return(XPtr<funcPtr>(new funcPtr(&rproposal_distr)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}
 */

/*
// [[Rcpp::export]]
NumericVector callViaXPtr(const NumericVector x, SEXP xpsexp, int n, double mult) {
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  NumericVector y = fun(x, n, mult);
  return (y);
}
*/
// [[Rcpp::export]]
List metropolisCbeta(NumericVector x, double beta, int num_iterations_mcmc,
                              SEXP rproposal_fn_xpsexp, SEXP dproposal_fn_xpsexp,
                           List other_params, int debug) {
  XPtr<rDistrFnPtr> xpfun_r(rproposal_fn_xpsexp);
  rDistrFnPtr rproposal_distr = *xpfun_r;
  
  XPtr<dCondDensityFnPtr> xpfun_d(dproposal_fn_xpsexp);
  dCondDensityFnPtr dproposal_distr = *xpfun_d;
  
  int n = x.size();
  //NumericVector proposal(n);
  NumericVector proposal;
  double log_old_p, log_new_p, log_diff;
  double log_old_g, log_new_g;
  double new_fa;
  bool do_compute_old_p=true;
  int lower_mult=10;
  int upper_mult=10000; 
  NumericMatrix all_x_vals;
  if(debug){

    NumericMatrix all_x_vals2(num_iterations_mcmc*4, n); //TODO 4 is hardcoded. change.
    all_x_vals= all_x_vals2;
  }
  
  int num_potential_changes=0;
  int num_actual_changes=0;
  int num_proposals_indomain = 0;
  int num_proposals = 0;
  //NumericVector covariances = NumericVector::create(0.01,0.05,0.1, 0.3);
  
  for(int i =1; i<=num_iterations_mcmc; i++){
    for(int j=lower_mult; j<=upper_mult; j=j*10){ // default: 10^3,10^5
      num_proposals++;
        /*NumericVector params1 = j*x;
        NumericVector params2 =j*(1-x);
        for(int index=0;index<n;index++){
          proposal[index] = rbeta(1, params1[index],params2[index])[0];
        }
        */
 
      //proposal = rproposal_distr(x, n, j);
      //XPtr<funcPtr> xpfun = XPtr<funcPtr>(new funcPtr(&rproposal_distr));//putFunPtrInXPtr(funname);
      //funcPtr fun = *xpfun;
      
      proposal = rproposal_distr(x,n,j);
      //proposal = callViaXPtr(x, proposal_fn,n,j);
      //std::cout<<x<<"ggf"<<proposal<<std::endl;
      //std::cout<<pow(exp(faC(x)),(1-beta));
      new_fa = faC(proposal);
      if(new_fa != -INFINITY){
        num_proposals_indomain ++;
        if(do_compute_old_p){
          log_old_p = faC(x)*(1-beta) + fbC(x, other_params)*beta;
          do_compute_old_p = false;
        }
        log_new_p = new_fa*(1-beta) + fbC(proposal, other_params)*beta;
        log_new_g = dproposal_distr(proposal, x, n, j);
        log_old_g = dproposal_distr(x, proposal, n, j);
        log_diff = (log_new_p - log_old_p) + (log_old_g -log_new_g);
        //std::cout<<new_p/old_p<<std::endl; 
        if(log_new_p != -INFINITY){
          num_potential_changes++;
          if(exp(log_diff) > runif(1)[0]){
            for(int index=0; index<n; index++){
              x[index] = proposal[index];
            }
            do_compute_old_p=true;
            num_actual_changes++;
            
          }
        } /*else {
          std::cout<<log_new_p<<new_fa<<" "<<fbC(proposal, other_params)<<log_old_p<<"  gf  "<<proposal<<std::endl;
        }*/
      }
      if(debug){
        all_x_vals(num_proposals-1,_) = x;
      }
    }
  }
  //std::cout<<num_potential_changes<<" Number of potential changes"<<std::endl;
  return List::create(
    _["x_vec"]=x,
    _["num_proposals"]=num_proposals,
    _["num_proposals_indomain"]=num_proposals_indomain,
    _["num_potential_changes"]=num_potential_changes,
    _["num_actual_changes"]=num_actual_changes,
    _["all_x_mat"] = all_x_vals
  );
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