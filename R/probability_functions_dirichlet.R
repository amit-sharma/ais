# Log-probability of "easy" distribution
# Uniform distribution
fa = function(e){
  if(any(e<0) || any(e>1))
    return(log(0))
  
  out <- 0 # uniform distribution log likelihood
  #out <- dhyperdirichlet(p,parameters, log=TRUE)
  #return(Jacobian(e)*exp(out))
  #return(log(Jacobian(e))+ out)
  return(out)
}

# Log-probability of "hard" distribution
#TODO making normalizing constant log scale to reduce overflow
log_likelihood <- function(theta, powers_dirichlet) {
  log_like = 0
  if(DISTR=="dirichlet"){
    log_like = sum(powers_dirichlet * log(theta))
    normalizing_constant = gamma(sum(powers_dirichlet+1))/prod(gamma(powers_dirichlet+1))
  } else if(DISTR=="hyperdirichlet") {
    theta_mat = matrix(theta, ncol=4,byrow=TRUE)
    theta_sum_terms = rowSums(theta_mat)
    log_like=sum(powers_dirichlet * log(theta_sum_terms))
    normalizing_constant = (gamma(sum(powers_dirichlet+4)) * gamma(4)^8)/prod(gamma(powers_dirichlet+4))
  }
  
  log_like = log_like + log(normalizing_constant)
  return(log_like)
} 

#TODO make sum of p equal to 1 as the condition by using float_equal
fb <- function(e, other_params){
  e <- c(1,e)
  p <- e_to_p(e)
  if(any(p<0) || any(p>1) || !isTRUE(all.equal.numeric(sum(p),1)))
    return(log(0))
  
  out <- log_likelihood(p, powers_dirichlet = other_params)
  #out <- dhyperdirichlet(p,parameters, log=TRUE)
  #return(Jacobian(e)*exp(out))
  print(paste(log(Jacobian(e)), out))
  return(log(Jacobian(e))+ out)
}

