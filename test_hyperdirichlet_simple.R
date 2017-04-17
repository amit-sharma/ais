# This example uses hyperdirichlet distribution.
# first testing 10 dimensional dirichlet distribution.

library(assertthat)
library(parallel)

library(Rcpp)
#library(inline)

library(hyperdirichlet)
library(ais)
library(SimplicialCubature)


set.seed(1)       # Seed the random number generator for reproducibility
DISTR="hyperdirichlet"
USE_CPP=TRUE
CUSTOM_PROPOSAL=FALSE
DO_PARALLEL=T

#dirichlet converges to exact at 10k and 200
K = 10000     # Number of annealed transitions per run, default 10000
replicates = 200  # Number of AIS runs, default 200


#TODO: check how approximation worsens off as degree increases.
Main <- function(){
  
  # Inverse temperatures
  betas = cooling2(K, exponent = -8) # seq(1, K)/K #
  
  if(DISTR=='dirichlet'){
    powers_dirichlet = c(1,rep(0,7))
    n=length(powers_dirichlet)-1 #            # Number of dimensions of distribution (1 less than 10)
    theta_sum_vec = list(
      "000"=c(0),
      "001"= c(1),
      "010"= c(2),
      "011"= c(3),
      "100"= c(4),
      "101"= c(5),
      "110"= c(6),
      "111"= c(7)
    )
    a=1
    b=1
  } else if (DISTR=="hyperdirichlet"){
    n=length(powers_dirichlet)*4 -1
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
    a=4
    b=4
  } else {
    print("Error in distribution")
  }
  
    
  
  # List of samples from the "easy" distribution
  # This is in the e domain (e_to_p converts it to probability domain)
  # so do not need to worry about summing to 1.
  samples = replicate(replicates, runif(n), simplify = FALSE)
  
  
  
  names(powers_dirichlet) = names(theta_sum_vec)
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
      transition = metropolisC2, 
      num_iterations_mcmc=10,
      proposal_sample_fn = rproposal,
      proposal_cond_density_fn = dproposal_cond, 
      added_e_power=0,
      other_params=list(powers=powers_dirichlet, param_structure=theta_sum_vec),
      parallel=DO_PARALLEL,
      num_cores=6
    )
    end_time=Sys.time()
    print(end_time-start_time)
    
    true_val = generalized_dirichlet_integral(powers_dirichlet, n+1, a,b, uselogscale=TRUE, added_constant=0, do_print=TRUE)
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
  
  #For comparison, we also look at integral method. 
  calc_hyperdirichlet_integral(x=NA, N_vec=powers_dirichlet, theta_vec=theta_sum_vec, all_theta_colindex=seq(0, n),  maxEval=50000, added_constant=0, do_cuhre=FALSE, do_simplex=TRUE)
  
  # Results -----------------------------------------------------------------
  return(ais_weightsC)
  
} 


calc_hyperdirichlet_integral <- function(x, N_vec, theta_vec, all_theta_colindex,  maxEval=0, added_constant=0, do_cuhre=NULL, do_simplex=TRUE){
  #n <- dim(x)-1  # "r-1" because this is the dimension of the integrand, not the number of p's.
  n <- length(all_theta_colindex)-1
  
  column_names = rev(names(all_theta_colindex))
  counter=1
  curr_p = NA 
  log_likelihood <- function(params) {
    curr_p <<- params
    counter <<- counter+sum(params)
    if(counter %% 10000 ==0) print(counter)
    p_params = params
    
    log_like = 0
    if(any(p_params<0) || any(p_params>1) || sum(p_params)>1)
      return(log(0))
    
    #      theta = c(p_params, 1-sum(p_params))
    theta=p_params
    names(theta) = column_names # initialize it here instead of test.r
    #print(theta)
    
    for(i in 1:length(N_vec)){
      zxy_val_count = N_vec[i]
      valid_theta_names = theta_vec[[names(N_vec)[i]]]
      log_like = log_like + zxy_val_count * log(sum(theta[valid_theta_names]))
    }
    names(log_like)=NULL
    return(log_like)
  }   
  
  upper_limit = rep(1,n)
  #upper_limit[c(3,8,13)] = 0.0000001
  
  if(do_cuhre){
    f <- function(e){
      e <- c(1,e)
      p <- e_to_p(e)
      out <- log_likelihood(p)
      #out <- dhyperdirichlet(p,parameters, log=TRUE)
      out <- out + added_constant
      return(Jacobian(e)*exp(out))
    }
    out <- adaptIntegrate(f,lowerLimit=rep(0,n),upperLimit=upper_limit,tol=1e-3, maxEval=maxEval)
    print(paste("Done 1 approx integral", out$integral))
    
    out2 <- cuhre(n,1, f,lower=rep(0,n),upper=upper_limit, flags=list(final=0), max.eval=1e5)
    out3 <- vegas(n, 1, f, lower=rep(0,n),upper=upper_limit, flags=list(final=0), max.eval=1e5)
    print(paste("Comparing adaptIntegrate and Cuhre", out$integral, out2$value, out3$value))
  }
  if(do_simplex){
    f2 <- function(p){
      #HD=x
      #parameters <- params(HD)
      
      p = c(p, 1-sum(p))
      out <- log_likelihood(p)
      #out <- dhyperdirichlet(p,parameters, log=TRUE)
      out <- out + added_constant
      return(exp(out))
    }
    out4 <- adaptIntegrateSimplex(f2, CanonicalSimplex(n), maxEvals=maxEval, tol=1e-10)
    print(out4$integral)
    #print(paste("Comparing adaptIntegrate and simplex", out$integral, out4$integral))
    print(paste("Simplex integral output", out4$integral))
  }
  print(paste("Added constant power of e (for hyperdirichlet):", added_constant))
  return(out4$integral)
}
