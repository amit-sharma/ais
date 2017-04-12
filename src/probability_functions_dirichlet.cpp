#include <cmath>
#include "probability_functions.h"

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
double logJacobianC(NumericVector e){
  double logres=0;
  int len_vec = e.size();
  for(int i=0; i<len_vec-1; i++){
    logres = logres +  (len_vec-i-1)*log(1-e[i]);
  }
  return(logres);
}


// [[Rcpp::export]]
double log_likelihoodC(NumericVector theta, IntegerVector powers_dirichlet, Rcpp::List theta_sum_list){
  double log_like, log_normalizing_constant=0, temp_sum=0, curr_num_vars=0;
  
  /*
  log_like = sum(as<NumericVector>(powers_dirichlet) * log(theta));
  log_normalizing_constant = lgamma(sum(powers_dirichlet+1)) - sum(lgamma(powers_dirichlet+1));
  */
  
  
  // For hyperdirichlet
  
  int num_terms = theta_sum_list.size();
  NumericVector theta_sum_terms(num_terms);
  int ii;
  for(ii=0; ii<num_terms;ii++){
    theta_sum_terms[ii] = 0;
    Rcpp::IntegerVector theta_indices = as<IntegerVector>(theta_sum_list[ii]);
    curr_num_vars = theta_indices.size();
    
    //For normalizing constant
    temp_sum += (powers_dirichlet[ii] + curr_num_vars );
    log_normalizing_constant +=  (lgamma(curr_num_vars) - lgamma(powers_dirichlet[ii]+curr_num_vars)); 
    
    for(int index=0;index < curr_num_vars;index++){
      theta_sum_terms[ii] += theta[theta_indices[index]];
    }
   }
  log_normalizing_constant += lgamma(temp_sum);
  
   /*
  for(int index=0;index < theta.size();index= index + 4){
    theta_sum_terms[index/4] = theta[index] + theta[index+1] + theta[index+2] + theta[index+3];
  }
  */
  log_like=sum(as<NumericVector>(powers_dirichlet) * log(theta_sum_terms));
  //NumericVector temp=log(theta_sum_terms);std::cout<<temp<<std::endl;
  
  //prior code
  /*log_normalizing_constant = lgamma(sum(powers_dirichlet+4)) + theta_sum_terms.size()*lgamma(4) - sum(lgamma(powers_dirichlet+4));*/ 
/*  for(ii=0;ii<num_terms;ii++){
    curr_num_vars = (as<IntegerVector>(theta_sum_list[ii])).size();
    temp_sum += (powers_dirichlet[ii] + curr_num_vars );
    log_normalizing_constant +=  (lgamma(curr_num_vars) - lgamma(powers_dirichlet[ii]+curr_num_vars)); 
  }
  log_normalizing_constant += lgamma(temp_sum);
 */
  
  if(!(log_like >-90)){
    //std::cout<<log_like<<powers_dirichlet<<std::endl;
    //std::cout<<theta_sum_terms<<theta<<std::endl;
  }
  
  // Commenting this out when we no longer need normalization..that is only needed for test_hyperdirichlet
  //log_like = log_like + log_normalizing_constant;
  return log_like;
}

//' @title
//' e_to_pC
//' @descriptinocp
//' Convert e to p representation.
//' 
//' @param e a vector
//' 
//' @details
//' \code{e_to_pC} takes a vector.
//' @export
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


//' @title
//' fbC
//' @descriptinocp
//' 
//' 
//' @param e a vector containing values of the parameters.
//' @param other_params a list containing other parameters
//' 
//' @details
//' \code{e_to_pC} takes a vector.
//' @export
// [[Rcpp::export]]
double fbC(NumericVector e, List other_params){
  //Already checked in faC
 if(is_true(any(e<0)))
    return -INFINITY ;
 if(is_true(any(e>1)))
      return -INFINITY;
 //if(is_true(any(e==1)))std::cout<<e<<std::endl;
  
  //processing list
  IntegerVector powers_dirichlet= as<IntegerVector>(other_params["powers"]);
  Rcpp:List theta_sum_list = other_params["param_structure"];
  
  double log_jacobian, out;//jacobian
   // because we know that powers will be integers
  NumericVector p;
  p = e_to_pC(e);
  
  if(is_true(any(p<0)))
    std::cout<<"ttr"<<e<<p<<std::endl;
  out = log_likelihoodC(p, powers_dirichlet, theta_sum_list);
  // Note: Jacobian will always be zero if one of the e is exactly 1.
  //jacobian = JacobianC(e);
  log_jacobian = logJacobianC(e);
  //if(jacobian==0) std::cout<<log_jacobian+out<<std::endl;
  //std::cout<<powers_dirichlet<<std::endl;
  //std::cout<<log(jacobian)<<" "<<out<<std::endl;
  return log_jacobian + out;
}
