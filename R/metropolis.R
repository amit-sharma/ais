metropolis = function(x, fa, fb, beta, jump, num_iterations_mcmc,  ...){
  for(i in 1:num_iterations_mcmc){
    for(j in c(1,3,10)){
      proposal = x + jump(j)
      #proposal = sapply(x, jump)
      old_p = exp(fa(x))^(1-beta) * exp(fb(x))^beta
      new_p = exp(fa(proposal))^(1-beta) * exp(fb(proposal))^beta
      
      if(new_p / old_p > runif(1)){
        x=proposal
      }
    }
  }
  return(x)
}
