# This example uses hyperdirichlet distribution.
# first testing 10 dimensional dirichlet distribution.

library(assertthat)
library(parallel)

library(hyperdirichlet)

source("R/ais.R")
source("R/metropolis.R")
source("R/log_pstar.R")
source("R/cooling.R")

set.seed(1)       # Seed the random number generator for reproducibility
DISTR="dirichlet"

powers_dirichlet= c(1, 2, 40, 5, 4, 12, 5, 11)
K = 10000        # Number of annealed transitions per run
replicates = 100  # Number of AIS runs

if(DISTR=='dirichlet'){
  n=length(powers_dirichlet)-1 #            # Number of dimensions of distribution (1 less than 10)
} else if (DISTR=="hyperdirichlet"){
  n=length(powers_dirichlet)*4 -1
} else {
  print("Error in distribution")
}

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
  #out <- dhyperdirichlet(p,parameters, log=TRUE)
  #return(Jacobian(e)*exp(out))
  return(log(Jacobian(e))+ out)
}



Main <- function(){
  # Inverse temperatures
  betas = cooling(K, exponent = -8)
  
  # List of samples from the "easy" distribution
  samples = replicate(replicates, runif(n), simplify = FALSE)
  
  # Collect importance weights using annealed importance sampling
  ais_weights = AIS(
    samples = samples, 
    betas = betas, 
    fa = fa, 
    fb = fb, 
    transition = metropolis, 
    #jump = function(){rt(n, (df + dfa) / 2)}
    jump = function(x){rnorm(n, mean=0, sd=0.05*x)},
    #jump = function(m){rbeta(1, m/(1-m), 1)},
    num_iterations_mcmc=10,
    parallel=FALSE
  )
  start_time=Sys.time()
  ais_weightsC = AISC(
    samples = samples, 
    betas = betas, 
    fa=faC,
    fb = fbC, 
    transition = metropolisC2, 
    num_iterations_mcmc=10,
    parallel=FALSE
  )
  end_time=Sys.time()
  print(end_time-start_time)
  
  # Collect importance weights using non-annealed importance sampling.
  # Note that this sampler collects K times more samples than the AIS sampler,
  # which should make them roughly similar in terms of total computational complexity.
  # The samples from AIS tend to be higher quality, even though they're each more 
  # difficult to compute.
  #s = matrix(rt(K * n * replicates, df = dfa), ncol = n)
  #is_weights = exp(rowSums(dt(s, df = df, log = TRUE) - dt(s, df = dfa, log = TRUE)))
  
  
  # Results -----------------------------------------------------------------
  
  print(paste("Integration (Mean of AIS weights)", mean(ais_weights)))  # Should be close to 1
  print(paste("Std Error of the mean", sd(ais_weights) / sqrt(length(ais_weights)))) # Estimated standard error (hopefully reliable)
  ## Adjusted sample size
  print(paste("Actual, Adjusted sample size", length(ais_weights), length(ais_weights)/(1+var(ais_weights))))
  
  print(paste("Integration (Mean of AIS weights)", mean(ais_weightsC)))  # Should be close to 1
  print(paste("Std Error of the mean", sd(ais_weightsC) / sqrt(length(ais_weightsC)))) # Estimated standard error (hopefully reliable)
  ## Adjusted sample size
  print(paste("Actual, Adjusted sample size", length(ais_weightsC), length(ais_weightsC)/(1+var(ais_weightsC))))
} 