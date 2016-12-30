cooling = function(K, exponent = -8){
  assert_that(exponent < 0)
  assert_that(K == as.integer(K))
  assert_that(K > 0)
  
  out = 1 - .95^(1:K / K * (exponent * log(10) / log(.95)))
  out[K] = 1
  
  out
}

cooling2 = function(K, exponent = -8){
  assert_that(exponent < 0)
  assert_that(K == as.integer(K))
  assert_that(K > 0)
  
  max_val = 0.999999
  min_val = 0.000001
  r = (max_val/min_val)^(1/K)
  out = min_val*(r^(1:K))
  out[K] = 1
  
  out
}