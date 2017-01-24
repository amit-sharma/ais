# This example uses hyperdirichlet distribution.
# first testing 10 dimensional dirichlet distribution.

library(assertthat)
library(parallel)

library(Rcpp)
library(inline)

library(hyperdirichlet)
library(ais)


set.seed(1)       # Seed the random number generator for reproducibility
DISTR="hyperdirichlet"
USE_CPP=TRUE
CUSTOM_PROPOSAL=FALSE
powers_dirichlet = c(1, 20, 400, 5, 60, 1,5000, 80)

K = 1000     # Number of annealed transitions per run, default 10000
replicates = 20  # Number of AIS runs, default 200

if(DISTR=='dirichlet'){
  n=length(powers_dirichlet)-1 #            # Number of dimensions of distribution (1 less than 10)
} else if (DISTR=="hyperdirichlet"){
  n=length(powers_dirichlet)*4 -1
} else {
  print("Error in distribution")
}


Main <- function(){
  # Inverse temperatures
  betas = cooling2(K, exponent = -8)
  
  # List of samples from the "easy" distribution
  samples = replicate(replicates, runif(n), simplify = FALSE)
  
  theta_sum_vec = list(
    "000"=c(0,1,2,3),
    "001"= c(4,5,6,7),
    "010"= c(8,9,10,11),
    "011"= c(12,13,14,15),
    "100"= c(16,17,18,19),
    "101"= c(20,21,22,23),
    "110"= c(24,25,26,27),
    "111"= c(28,29,30,31)
  )
  if(!USE_CPP){
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
      num_iterations_mcmc=50,
      other_params=powers_dirichlet,
      parallel=FALSE
    )
    print(paste("Integration (Mean of AIS weights)", mean(ais_weights)))  # Should be close to 1
    print(paste("Std Error of the mean", sd(ais_weights) / sqrt(length(ais_weights)))) # Estimated standard error (hopefully reliable)
    ## Adjusted sample size
    print(paste("Actual, Adjusted sample size", length(ais_weights), length(ais_weights)/(1+var(ais_weights))))
    
  } else {
    start_time=Sys.time()
    
    if(CUSTOM_PROPOSAL==TRUE){
      sourceCpp("distr_functions.cpp")
      rproposal = putFunPtrInXPtrCustom()
    } else {
      rproposal = putFunPtrInXPtr("default")  #FROM AIS
      dproposal_cond = getCondDensityFuncXPtr("default")
    }
    
    ais_weightsC = AISC(
      samples = samples, 
      betas = betas, 
      fa=faC,
      fb = fbC, 
      transition = metropolisCbeta, 
      num_iterations_mcmc=10,
      proposal_sample_fn = rproposal,
      proposal_cond_density_fn = dproposal_cond, 
      other_params=list(powers=powers_dirichlet, param_structure=theta_sum_vec),
      parallel=F,
      num_cores=6
    )
    end_time=Sys.time()
    print(end_time-start_time)
  
    print(paste("Integration (Mean of AIS weights)", mean(ais_weightsC), mean(ais_weightsC, na.rm=TRUE)))  # Should be close to 1
    print(paste("Std Error of the mean", sd(ais_weightsC) / sqrt(length(ais_weightsC)), sd(ais_weightsC, na.rm=TRUE) / sqrt(length(ais_weightsC)) )) # Estimated standard error (hopefully reliable)
    ## Adjusted sample size
    print(paste("Actual, Adjusted sample size", length(ais_weightsC), length(ais_weightsC)/(1+var(ais_weightsC))))
  } 
  
  
  
  # Collect importance weights using non-annealed importance sampling.
  # Note that this sampler collects K times more samples than the AIS sampler,
  # which should make them roughly similar in terms of total computational complexity.
  # The samples from AIS tend to be higher quality, even though they're each more 
  # difficult to compute.
  #s = matrix(rt(K * n * replicates, df = dfa), ncol = n)
  #is_weights = exp(rowSums(dt(s, df = df, log = TRUE) - dt(s, df = dfa, log = TRUE)))
  
  
  # Results -----------------------------------------------------------------
  return(ais_weightsC)
  
} 
