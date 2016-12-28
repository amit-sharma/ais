# This example uses hyperdirichlet distribution.
# first testing 10 dimensional dirichlet distribution.

library(assertthat)
library(parallel)

library(hyperdirichlet)

set.seed(1)       # Seed the random number generator for reproducibility

n = 9            # Number of dimensions of distribution (1 less than 10)
powers_dirichlet= c(1, 2, 4, 5, 6, 1, 5, 8, 7,1)
K = 10000        # Number of annealed transitions per run
replicates = 100  # Number of AIS runs


# Log-probability of "easy" distribution
fa = function(e){
  e <- c(1,e)
  p <- e_to_p(e)
  if(any(p<0) || any(p>1) || sum(p)>1)
    return(log(0))
  
  out <- 0 # uniform distribution log likelihood
  #out <- dhyperdirichlet(p,parameters, log=TRUE)
  #return(Jacobian(e)*exp(out))
  #return(log(Jacobian(e))+ out)
  return(out)
}

# Log-probability of "hard" distribution
log_likelihood <- function(params) {
  p_params = params
  
  log_like = 0
  
  theta=p_params
  log_like = powers_dirichlet %*% log(theta)
  normalizing_constant = gamma(sum(powers_dirichlet+1))/prod(gamma(powers_dirichlet+1))
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
  num_iterations_mcmc=10
  #parallel=FALSE
)

# Collect importance weights using non-annealed importance sampling.
# Note that this sampler collects K times more samples than the AIS sampler,
# which should make them roughly similar in terms of total computational complexity.
# The samples from AIS tend to be higher quality, even though they're each more 
# difficult to compute.
#s = matrix(rt(K * n * replicates, df = dfa), ncol = n)
#is_weights = exp(rowSums(dt(s, df = df, log = TRUE) - dt(s, df = dfa, log = TRUE)))


# Results -----------------------------------------------------------------

mean(ais_weights)  # Should be close to 1
#mean(is_weights)   # Should also be close to 1, but will usually be too low
########
sd(ais_weights) / sqrt(length(ais_weights)) # Estimated standard error (hopefully reliable)
#sd(is_weights) / sqrt(length(is_weights))   # Estimated standard error (unreliable)

## Adjusted sample size
length(ais_weights)/(1+var(ais_weights))
