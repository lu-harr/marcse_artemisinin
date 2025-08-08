# pinching from idem-lab/ir_cube:

# given a vector of sample sizes N, predicted proportions p, and overdispersion
# parameter rho (on unit interval, where 0 means no overdispersion and 1 means
# maximum), returna betatbinomial-distributed variable greta array

betabinomial_p_rho <- function(N, p, rho) {
  
  # model the observation (betabinomial) sd as a multiplier on the binomial sd,
  # accounting for additional error due to nonindependent sampling of individuals
  # from the population. This is based on the INLA parameterisation
  
  # solve for a and b:
  #   p = a / (a + b)
  #   rho = 1 / (a + b + 1)
  a <- p * (1 / rho - 1)
  b <- a * (1 - p) / p
  
  # define betabinomial according to the greta interface
  greta::beta_binomial(size = N, alpha = a, beta = b)
  
}
