library(greta.gp) # we won't need my diy version for this one
source("code/build_design_matrix.R")

# a start - works okay
coords <- mut_data %>%
  dplyr::select(x, y, year, scaled_year)

define_greta_parameters <- function(){
  list(
    beta = normal(0, 3, dim=3), # I have three betas so dim is 3 (intercept, time, PfPR)
    kernel_lengthscale_space = normal(0, 3, truncation = c(0, Inf)),
    kernel_lengthscale_time = normal(0, 3, truncation = c(0, Inf)),
    kernel_sd = normal(0, 1, truncation = c(0, Inf))
  )
}

parameters <- define_greta_parameters() # are there any choices that should be arguments?


kernel <- mat52(lengthscales = c(parameters$kernel_lengthscale_space, 
                                 parameters$kernel_lengthscale_space, 
                                 parameters$kernel_lengthscale_time),
                variance = parameters$kernel_sd ^ 2)

X_obs <- build_design_matrix(covariates, 
                             coords, 
                             scale = FALSE, 
                             temporal_range = pfpr_years)
# we want this to give us back scaled years but need to provide unscaled years
# so that it can grab from correct annual covt raster
message(nrow(X_obs))

coords <- coords %>%
  dplyr::select(x, y, scaled_year)

random_field <- gp(x = coords, # lon, lat, year
                   kernel = kernel,
                   inducing = unique(coords)  # shouldn't do much to simulated data but might be very helpful for actual data
)

gp_mean_obs <- X_obs %*% parameters$beta + random_field

mut_prob_obs <- ilogit(gp_mean_obs)

distribution(mut_data$present) <- greta::binomial(size = mut_data$tested,
                                                  prob = mut_prob_obs)

beta = parameters$beta
m <- model(parameters$kernel_lengthscale_space, 
           parameters$kernel_lengthscale_time, 
           parameters$kernel_sd, 
           parameters$beta)

draws <- mcmc(m, n_samples = 3000, 
              initial_values = replicate(4, initials(beta = rep(0,3)), 
                                         simplify = FALSE))

png("output/mat52_trace.png", height=750, width=1500)
bayesplot::mcmc_trace(draws)
dev.off()

r_hats <- coda::gelman.diag(draws,
                            autoburnin = FALSE,
                            multivariate = FALSE)
summary(r_hats$psrf)

# save everything and do the prediction in a separate script
write_rds(mut_data, "output/mat52_model/mut_data.rds")
write_rds(parameters, "output/mat52_model/parameters.rds")
# write_rds(mut_prob_obs, "output/mat52_model/mut_prob_obs.rds")
write_rds(kernel, "output/mat52_model/kernel.rds")
write_rds(random_field, "output/mat52_model/random_field.rds")
write_rds(m, "output/mat52_model/m.rds")
write_rds(draws, "output/mat52_model/draws.rds")









