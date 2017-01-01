AIS = function(samples, betas, fa, fb, transition, num_iterations_mcmc,parallel = TRUE, num_cores=1, other_params=NULL,  ...){
  run_fn <- function(x){
    run(x, betas = betas, fa = fa, fb = fb, transition = transition, num_iterations_mcmc=num_iterations_mcmc, 
        other_params=other_params, ...)
  }
  if(parallel){
    res = mclapply(samples,run_fn, mc.cores=num_cores)
  }else{
    res = lapply(samples, run_fn)
  }
  
  simplify2array(res)
}

run = function(x, betas, fa, fb, transition, other_params,  ...){
  K = length(betas)
  
  assert_that(all(betas <= 1))
  assert_that(all(betas >= 0))
  assert_that(all(betas == sort(betas)))
  assert_that(betas[K] == 1)
  
  # Empty vectors for storing each sample's negative energy under 
  # both parent distributions. Throughout, "a" refers to the prior/simple distribution
  # and "b" refers to the intractable one.
  f_as = numeric(K)
  f_bs = numeric(K)
  
  for(k in 1:K){
    # Sample at new temperature
    x = transition(x, fa, fb, betas[k], ...)
    
    # save negative energies under both distributions
    f_as[k] = fa(x)
    f_bs[k] = fb(x,other_params)
  }
  
  # Betas in numerator goes from 1:K
  # Betas in denominator go from 0:(K-1)
  w = exp(
    sum(
      log_pstar(f_as, f_bs, betas) - log_pstar(f_as, f_bs, c(0, betas[-K]))
    )
  )
  #print("Completed run")
  w
}



AISC <- function(samples, betas, fa, fb, transition, num_iterations_mcmc, other_params=NULL, 
                parallel = TRUE, num_cores=1, ...){
  run_fn <- function(x){
    runC(x, betas = betas, fa = fa, fb = fb, transition = transition,num_iterations_mcmc=num_iterations_mcmc,
         other_params=other_params, ...)
  }
  if(parallel){
    res = mclapply(samples,run_fn, mc.cores=num_cores)
  }else{
    res = lapply(samples, run_fn)
  }
  
  simplify2array(res)
  # if(parallel){
  #   f = mclapply
  # }else{
  #   f = lapply
  # }
  # 
  # simplify2array(
  #   f(
  #     samples,
  #     function(x){
  #       runC(x, betas = betas, fa = fa, fb = fb, transition = transition, num_iterations_mcmc=num_iterations_mcmc, ...)
  #     },
  #     mc.cores=4
  #   )
  # )
}

runC = function(x, betas, fa, fb, transition, num_iterations_mcmc, other_params, ...){
  K = length(betas)
  
  assert_that(all(betas <= 1))
  assert_that(all(betas >= 0))
  assert_that(all(betas == sort(betas)))
  assert_that(betas[K] == 1)
  
  # Empty vectors for storing each sample's negative energy under 
  # both parent distributions. Throughout, "a" refers to the prior/simple distribution
  # and "b" refers to the intractable one.
  f_as = numeric(K)
  f_bs = numeric(K)
  
  for(k in 1:K){
    # Sample at new temperature
    #x = transition(x, fa, fb, betas[k], ...)
    x = transition(x, betas[k], num_iterations_mcmc, other_params, ...)
    
    # save negative energies under both distributions
    f_as[k] = fa(x)
    f_bs[k] = fb(x,other_params)
  }
  
  # Betas in numerator goes from 1:K
  # Betas in denominator go from 0:(K-1)
  w = exp(
    sum(
      log_pstar(f_as, f_bs, betas) - log_pstar(f_as, f_bs, c(0, betas[-K]))
    )
  )
  #print("Completed run")
  w
}
