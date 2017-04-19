AIS = function(samples, betas, fa, fb, transition, num_iterations_mcmc,parallel = TRUE, num_cores=1, added_e_power=0, other_params=NULL,  ...){
  run_fn <- function(x){
    run(x, betas = betas, fa = fa, fb = fb, transition = transition, num_iterations_mcmc=num_iterations_mcmc, 
        added_e_power=added_e_power,
        other_params=other_params, ...)
  }
  if(parallel){
    res = mclapply(samples,run_fn, mc.cores=num_cores)
  }else{
    res = lapply(samples, run_fn)
  }
  
  simplify2array(res)
}

run = function(x, betas, fa, fb, transition, added_e_power, other_params,  ...){
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
    x = transition(x, fa, fb, betas[k], other_params, ...)
    
    # save negative energies under both distributions
    f_as[k] = fa(x)
    f_bs[k] = fb(x,other_params)
    print(paste(f_as[k],f_bs[k]))
  }
  
  # Betas in numerator goes from 1:K
  # Betas in denominator go from 0:(K-1)
  w = exp(
    sum(
      log_pstar(f_as, f_bs, betas) - log_pstar(f_as, f_bs, c(0, betas[-K]))
    ) + added_e_power
  )
  #print("Completed run")
  return(w)
}



#' Title
#'
#' @param samples 
#' @param betas  Temperature vector for annealing
#' @param fa 
#' @param fb 
#' @param transition 
#' @param num_iterations_mcmc 
#' @param proposal_sample_fn 
#' @param proposal_cond_density_fn 
#' @param other_params 
#' @param parallel 
#' @param num_cores 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
AISC <- function(samples, betas, fa, fb, transition, num_iterations_mcmc, 
                 proposal_sample_fn, proposal_cond_density_fn, 
                 added_e_power,
                 other_params=NULL, 
                parallel = TRUE, num_cores=1, ...){
  run_fn <- function(x){
    runC(x, betas = betas, fa = fa, fb = fb, transition = transition,
         added_e_power=added_e_power,
         num_iterations_mcmc=num_iterations_mcmc,
         proposal_sample_fn=proposal_sample_fn,
         proposal_cond_density_fn= proposal_cond_density_fn,
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

runC = function(x, betas, fa, fb, transition, 
                added_e_power,
                num_iterations_mcmc, 
                proposal_sample_fn, proposal_cond_density_fn,
                other_params, ...){
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
  
  # plot data matrix
  x_plot_matrix = matrix(ncol=length(x)+1, nrow=50000)
  plot_matrix_i = 1
  
  num_changes = 0
  #print(paste("initial x", x))
  sum_num_actual_changes = 0
  sum_num_proposals = 0
  sum_num_proposals_indomain=0
  sum_num_potential_changes=0
  for(k in 1:K){
    x_plot_matrix[plot_matrix_i,]=e_to_pC(x)
    plot_matrix_i = plot_matrix_i+1
    
    # Sample at new temperature
    #x = transition(x, fa, fb, betas[k], ...)
    # This function changes x in place as well
    ret_list = transition(x, betas[k], num_iterations_mcmc, 
                   proposal_sample_fn, proposal_cond_density_fn, other_params, ...)
    x = ret_list[["x_vec"]]
    sum_num_actual_changes = sum_num_actual_changes+ret_list[["num_actual_changes"]]
    sum_num_potential_changes = sum_num_potential_changes+ret_list[["num_potential_changes"]]
    sum_num_proposals = sum_num_proposals+ret_list[["num_proposals"]]
    sum_num_proposals_indomain = sum_num_proposals_indomain+ret_list[["num_proposals_indomain"]]
    
    
    #if(){ num_changes = num_changes+ 1}
    
    # save negative energies under both distributions
    f_as[k] = fa(x)
    f_bs[k] = fb(x,other_params)
    #print(paste(f_as[k],f_bs[k]))a
    
  }
  plot_matrix_fn(x_plot_matrix) 
  print(paste("Acceptance ratio", sum_num_actual_changes/sum_num_proposals, sum_num_potential_changes/sum_num_proposals, sum_num_proposals_indomain/sum_num_proposals, sum_num_proposals ))
  # Betas in numerator goes from 1:K
  # Betas in denominator go from 0:(K-1)
  w = exp(
    sum(
      log_pstar(f_as, f_bs, betas) - log_pstar(f_as, f_bs, c(0, betas[-K]))
    ) +
      added_e_power
  )
  #print("Completed run")
  return(w)
}


plot_matrix_fn <- function(mat){
  df = as.data.frame(mat)
  df = df %>% mutate(iter=row_number()) %>% filter(!is.na(V1))
  df2 = gather(df, param_name, param_value, starts_with("V"))
  p1=ggplot(df2, aes(x=iter, y=param_value)) + geom_line() + facet_wrap(~param_name)
  print(p1)
  p2 = ggplot(df2, aes(x=param_value)) + geom_density()+facet_wrap(~param_name)
  print(p2)
  df3 = group_by(df2, param_name) %>% mutate(running_mean = cumsum(param_value)/iter) %>%
    ungroup()
  p3 = ggplot(df3, aes(x=iter, y=running_mean)) + geom_line() + facet_wrap(~param_name)
  print(p3)
  
}