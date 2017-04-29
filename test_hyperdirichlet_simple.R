# This example uses hyperdirichlet distribution.
# first testing 10 dimensional dirichlet distribution.

# WARN: powers_dirichlet and theta_sum_vec should be same length and indices in the same order.

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

library(assertthat)
library(parallel)

library(Rcpp)
#library(inline)

library(hyperdirichlet)
library(ais)
library(SimplicialCubature)


set.seed(1)       # Seed the random number generator for reproducibility
#DISTR="hyperdirichlet"
USE_CPP=TRUE
CUSTOM_PROPOSAL=FALSE
DO_PARALLEL=T
DEBUG=F

#dirichlet converges to exact at 10k and 200
K = 20000   # Number of annealed transitions per run, default 10000
replicates = 200  # Number of AIS runs, default 200
num_powers=8 # 32 for dirichlet

#TODO: check how approximation worsens off as degree increases.
Main <- function(DISTR, NUM_CORES=6){
  
  # Inverse temperatures
  betas = cooling2(K, exponent = -8) # 
  #betas= seq(1, K)/K 
  
  if(DISTR=='dirichlet'){
    powers_dirichlet = round(runif(num_powers)*5)
    n=length(powers_dirichlet)-1 #            # Number of dimensions of distribution (1 less than 10)
    # theta_sum_vec_all = list(
    #   "000"=c(0),
    #   "001"= c(1),
    #   "010"= c(2),
    #   "011"= c(3),
    #   "100"= c(4),
    #   "101"= c(5),
    #   "110"= c(6),
    #   "111"= c(7),
    #   "1000"=c(8),
    #   "1001"=c(9),
    #   "1010"=c(10),
    #   "1011"=c(11),
    #   "1100"=c(12),
    #   "1101"=c(13),
    #   "1110"=c(14),
    #   "1111"=c(15)
    # )
    # theta_sum_vec=theta_sum_vec_all[1:num_powers]
    theta_sum_vec = lapply(seq(0,n), function(x) c(x))
    names(theta_sum_vec)=seq(0,n)
    a=1
    b=1
  } else if (DISTR=="hyperdirichlet"){
    powers_dirichlet = round(runif(num_powers)*50)
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
  } else if(DISTR=="hyper-data") {
    zxy_data = as.data.frame(matrix(
      c(
        rep(c(0,0,0), 2),
        rep(c(0,0,1), 10),
        rep(c(1, 1, 1), 1)
      ),
      byrow = TRUE,ncol=3)
    ) %>%
      rename(z_vec=V1, x_vec=V2, y_vec=V3)
    
    freq_unique_zxy = group_by(zxy_data, z_vec, x_vec, y_vec) %>% 
      summarize(num_data_points=n()) %>%
      ungroup() %>%
      mutate(rowid_str=paste0(z_vec, x_vec, y_vec))
    param_vec = freq_unique_zxy$num_data_points
    names(param_vec)=freq_unique_zxy$rowid_str
    

    
    r = 24#nrow(freq_unique_zxy)
    #prior = factorial(r-1)  # CHECK this ..what is r here, number of parameters+-1
    #log_prior=log(prior)
    a=4
    b=2
    # theta_vec = list(
    #   "000"=c(0,1,2,3),
    #   "001"= c(4,5,6,7),
    #   "010"= c(8,9),
    #   "011"= c(10,11),
    #   "100"= c(12,13),
    #   "101"= c(14,15),
    #   "110"= c(16,17,18,19),
    #   "111"= c(20,21,22,23)
    # )
    
    
    
    
    allowed_rx_vals = c(0, 2, 3)
    allowed_ry_vals = c(0, 1, 2, 3)
    # if specify a non-uniform prior, then this order will matter.
    all_unique_thetas_vec = sort(
      outer(
        c("0","1"),
        c(outer(allowed_rx_vals, allowed_ry_vals, paste, sep="_")),
        paste, sep="_"
      ) # actual order should not matter as long as we follow the same order throughout this function
    )
    
    #r = 2*length(allowed_rx_vals)*length(allowed_ry_vals)# number of parameters; number of free parameters=23 rx={0,2,3}x ry={0,1,2,3} #xrow(distinct(zxy_data, z_vec, x_vec, y_vec)
    
    theta_y_vec = list("000"=c(0,2),
                       "001"=c(1,3),
                       "010"=c(0,1),
                       "011"=c(2,3),
                       "100"=c(0,2),
                       "101"=c(1,3),
                       "110"=c(0,1),
                       "111"=c(2,3)
    )
    theta_vec = list(
      "000" = paste("0", c( t(outer(c(0,2), theta_y_vec[["000"]], paste, sep="_"))), sep="_"),
      "001" = paste("0", c( t(outer(c(0,2), theta_y_vec[["001"]], paste, sep="_")) ), sep="_"),
      "010" = paste("0", c( t(outer(c(3), theta_y_vec[["010"]], paste, sep="_")) ), sep="_"),
      "011" = paste("0", c( t(outer(c(3), theta_y_vec[["011"]], paste, sep="_")) ), sep="_"),
      "100" = paste("1", c( t(outer(c(0), theta_y_vec[["100"]], paste, sep="_")) ), sep="_"),
      "101" = paste("1", c( t(outer(c(0), theta_y_vec[["101"]], paste, sep="_")) ), sep="_"),
      "110" = paste("1", c( t(outer(c(2,3), theta_y_vec[["110"]], paste, sep="_")) ), sep="_"),
      "111" = paste("1", c( t(outer(c(2,3), theta_y_vec[["111"]], paste, sep="_")) ), sep="_")
    )
    theta_vec_ais = sapply(theta_vec, function(x){sapply(x, function(y){which(all_unique_thetas_vec==y)-1})}) # converting to indices, -1 because c++ indices start at 0
    
    #anum_vec = c(4, 4, 2, 2, 2, 2, 4, 4)
    anum_vec = sapply(theta_vec, length)

    
    param_vec2 = rep(0, length(theta_vec_ais))
    names(param_vec2)=names(theta_sum_vec)
    for(index_id in names(param_vec)){
      param_vec2[index_id]= param_vec[index_id]
    }
    powers_dirichlet = param_vec2
    n=length(powers_dirichlet)*a/2 +length(powers_dirichlet)*b/2 -1
    theta_sum_vec = theta_vec_ais
   # a=4
  #  b=2
    
  } else {
    print("Error in distribution")
  }
  
    
  
  # List of samples from the "easy" distribution
  # This is in the e domain (e_to_p converts it to probability domain)
  # so do not need to worry about summing to 1.
  #samples = replicate(replicates, runif(n), simplify = FALSE)
  orig_samples = replicate(replicates, runif(n+1), simplify = FALSE)
  p_samples=lapply(orig_samples, function(elem_arr){ elem_arr/sum(elem_arr)})
  samples = lapply(p_samples, function(elem_arr) {p_to_e(elem_arr)[-1]})
  
 
  added_constant= sum(powers_dirichlet) 
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
      transition = metropolisCbeta, 
      num_iterations_mcmc=20, # was 40.
      proposal_sample_fn = rproposal,
      proposal_cond_density_fn = dproposal_cond, 
      added_e_power=added_constant,
      other_params=list(powers=powers_dirichlet, param_structure=theta_sum_vec),
      parallel=DO_PARALLEL,
      num_cores=NUM_CORES,
      debug=DEBUG
    )
    end_time=Sys.time()
    print(end_time-start_time)
    
    #true_val = generalized_dirichlet_integral(powers_dirichlet, n+1, a,b, uselogscale=TRUE, added_constant=0, do_print=TRUE)
    print(paste("ADDED CONSTANT", added_constant))
    true_val_byformula= dirichlet_integral_formula(DISTR, powers_dirichlet, theta_sum_vec, added_constant=added_constant) 
    val_by_ais = mean(ais_weightsC)
    print(paste("True exact value", true_val_byformula))
    print(paste("Integration (Mean of AIS weights)", val_by_ais, mean(ais_weightsC, na.rm=TRUE), log(mean(ais_weightsC))))  # Should be close to 1
    print(paste("Std Error of the mean", sd(ais_weightsC) / sqrt(length(ais_weightsC)), sd(ais_weightsC, na.rm=TRUE) / sqrt(length(ais_weightsC)) )) # Estimated standard error (hopefully reliable)
    ## Adjusted sample size
    print(paste("Actual, Adjusted sample size", length(ais_weightsC), length(ais_weightsC)/(1+var(ais_weightsC/val_by_ais))))
    print(paste("Percent error", mean(ais_weightsC)/true_val_byformula*100 - 100))
  } 
  
  
  
  # Collect importance weights using non-annealed importance sampling.
  # Note that this sampler collects K times more samples than the AIS sampler,
  # which should make them roughly similar in terms of total computational complexity.
  # The samples from AIS tend to be higher quality, even though they're each more 
  # difficult to compute.
  #s = matrix(rt(K * n * replicates, df = dfa), ncol = n)
  #is_weights = exp(rowSums(dt(s, df = df, log = TRUE) - dt(s, df = dfa, log = TRUE)))
  
  #For comparison, we also look at integral method. 
  all_theta_colindex=seq(0,n);names(all_theta_colindex)=seq(0,n)
  calc_hyperdirichlet_integral(x=NA, N_vec=powers_dirichlet, theta_vec=theta_sum_vec, all_theta_colindex=all_theta_colindex,  maxEval=50000, added_constant=added_constant, do_cuhre=FALSE, do_simplex=TRUE)
  
  # Results -----------------------------------------------------------------
  return(ais_weightsC)
  
} 


dirichlet_integral_formula <- function(DISTR, powers_dirichlet, theta_sum_vec, added_constant){
  if(DISTR=="dirichlet"){
    log_ans = sum(lgamma(powers_dirichlet+1)) - lgamma(length(powers_dirichlet)+sum(powers_dirichlet)) + added_constant
  } else if (DISTR=="hyperdirichlet"){
    log_ans = sum(lgamma(powers_dirichlet+4) - lgamma(4)) - lgamma(4*length(powers_dirichlet)+sum(powers_dirichlet)) + added_constant
  } else if (DISTR=="hyper-data"){
    powers1= powers_dirichlet[c("000", "001", "110", "111")]
    powers2= powers_dirichlet[c("010", "011", "100", "101")]
    log_ans = sum(lgamma(powers1+4) - lgamma(4)) + sum(lgamma(powers2+2) - lgamma(2)) - lgamma(length(unlist(theta_sum_vec))+sum(powers_dirichlet)) + added_constant 
  }
  return(exp(log_ans))
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
