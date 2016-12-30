library(ais)
library(MCMCpack)
library(hyperdirichlet)

library(testthat)

#context("Evaluate probability functions proportional to density")
VEC_SIZE=10
DISTR="dirichlet"


e_prob_vec = runif(runif(1, min=1, max=VEC_SIZE))
e_neg_vec = runif(runif(1, min=1, max=VEC_SIZE), min=-0.1, max=0)
e_pos_vec = runif(runif(1, min=1, max=VEC_SIZE), min=1, max=1.1)

powers_dirichlet_e = sample(1:15, length(e_prob_vec)+1)


p_prob_vec = runif(runif(1, min=1, max=VEC_SIZE))
p_simplex_vec = p_prob_vec/sum(p_prob_vec)
p_neg_vec = runif(runif(1, min=1, max=VEC_SIZE), min=-0.1, max=0)
p_pos_vec = runif(runif(1, min=1, max=VEC_SIZE), min=1, max=1.1)

powers_dirichlet_p= sample(1:15, length(p_simplex_vec))

test_that("Rcpp prior is uniform",{
  expect_equal(fa(e_prob_vec), 0)
  expect_equal(fa(e_neg_vec), -Inf)
  expect_equal(fa(e_pos_vec), -Inf)
  expect_equal(fa(e_prob_vec), faC(e_prob_vec))
  expect_equal(fa(e_neg_vec), faC(e_neg_vec))
  expect_equal(fa(e_pos_vec), faC(e_pos_vec))
})

# sometimes the test may fail if ddirichlet from MCMCpack overflows. Be careful.
test_that("Rcpp likelihood calculation is correct",{
  expect_equal(log_likelihood(p_simplex_vec, powers_dirichlet_p), 
               log(ddirichlet(p_simplex_vec, powers_dirichlet_p+1)))
  expect_equal(log_likelihood(p_simplex_vec, powers_dirichlet_p), 
               log_likelihoodC(p_simplex_vec, powers_dirichlet_p))
})

test_that("Rcpp Jacobian is correct",{
  expect_equal(Jacobian(c(1,e_prob_vec)), JacobianC(e_prob_vec))
})

test_that("Rcpp e_to_p is correct",{
  expect_equal(e_to_p(c(1, e_prob_vec)), e_to_pC(e_prob_vec))
})

test_that("Rcpp likelihood (with Jacobian) is correct", {
  expect_equal(fb(e_prob_vec, powers_dirichlet_e), log(Jacobian(c(1, e_prob_vec)))+log_likelihood(e_to_p(c(1,e_prob_vec)),powers_dirichlet_e))
  expect_equal(fb(e_neg_vec, powers_dirichlet_e), -Inf)
  expect_equal(fb(e_pos_vec, powers_dirichlet_e), -Inf)
  expect_equal(fb(e_prob_vec, powers_dirichlet_e), fbC(e_prob_vec, powers_dirichlet_e))
  expect_equal(fbC(e_neg_vec, powers_dirichlet_e), -Inf)
  expect_equal(fbC(e_pos_vec, powers_dirichlet_e), -Inf)
})
