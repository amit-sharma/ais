metropolis = function(x, fa, fb, beta, jump, num_iterations_mcmc,  ...){
  do_change_old_p=TRUE
  log_old_p=NA
  log_new_p=NA
  for(i in 1:num_iterations_mcmc){
    for(j in c(1,3,9)){
      proposal = x + jump(j)
      #print(x)
      #print(proposal)
      #proposal = sapply(x, jump)
      log_new_fa =fa(proposal)
      if(is.finite(log_new_fa)){
        if(do_change_old_p){
          log_old_p = fa(x)*(1-beta) + fb(x)*beta
          do_change_old_p=FALSE
        }
        log_new_fb = fb(proposal)
        log_new_p = log_new_fa*(1-beta) + log_new_fb*beta
        log_diff = log_new_p - log_old_p      
        if(exp(log_diff) > runif(1)){
          x=proposal
          do_change_old_p = TRUE
        }
      }
    }
  }
  return(x)
}

