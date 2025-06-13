#setwd("~/Desktop/MARCSE/k13_seafrica")
source("code/setup.R")

# or install local version ...:
devtools::install_github("lu-harr/greta.gp.st.on.earth")
library(greta.gp)

# here are the bits I've added
?circmat
?degrees_to_radians
?great_circle_dist

# model with circular matern kernel over space?
coords <- mut_data %>%
  dplyr::select(x, y, year, scaled_year)

source("code/build_design_matrix.R")
X_obs <- build_design_matrix(covariates, 
                             coords, 
                             scale = FALSE, 
                             temporal_range = pfpr_years)
# we want this to give us back scaled years but need to provide unscaled years
# so that it can grab from correct annual covt raster

# ggplot() +
#   geom_sf(data = filter(world, continent == "Africa"), fill = "white") +
#   geom_point(data = coords, aes(x, y))
# 
# world_radians <- degrees_to_radians(world %>% filter(continent == "Africa") %>% st_geometry())

coords <- coords %>%
  # convert to radians - need to do this after covariate extraction above ..
  # check this doesn't mess things up when we get to prediction ...
  # will probably just need a radians_to_degrees?
  mutate(x = degrees_to_radians(x),
         y = degrees_to_radians(y))

# this is reassuring:
# ggplot() +
#   geom_sf(data = world_radians, fill = "white") +
#   geom_point(data = coords, aes(x, y))

# define_greta_parameters <- function(){
#   list(
#     beta = normal(0, 3, dim=3), # I have three betas so dim is 3 (intercept, time, PfPR)
#     kernel_lengthscale_space = normal(0, 0.1, truncation = c(0, Inf)),
#     kernel_lengthscale_time = normal(0, 3, truncation = c(0, Inf)),
#     kernel_sd = normal(0, 1, truncation = c(0, Inf))
#   )
# }
# 
# parameters <- define_greta_parameters() # are there any choices that should be arguments?

beta = normal(0, 3, dim=3) # I have three betas so dim is 3 (intercept, time, PfPR)
kernel_lengthscale_space = normal(0, 0.1, truncation = c(0, Inf))
kernel_lengthscale_time = normal(0, 3, truncation = c(0, Inf))
kernel_sd = normal(0, 1, truncation = c(0, Inf))

# kernel <- circmat(lengthscale = parameters$kernel_lengthscale_space, 
#                   variance = parameters$kernel_sd**2, 
#                   columns = c(1, 2),
#                   circumference = 1) + 
#           expo(lengthscales = parameters$kernel_lengthscale_time,
#                variance = parameters$kernel_sd**2,
#                columns = 3)

kernel <- circmat(lengthscale = kernel_lengthscale_space, 
                  variance = kernel_sd**2, 
                  columns = c(1, 2),
                  circumference = 1) + 
  expo(lengthscales = kernel_lengthscale_time,
       variance = kernel_sd**2,
       columns = 3)

# need to incorporate `kernel _sd` .. could go in both variance spots ..?

coords <- coords %>%
  dplyr::select(x, y, scaled_year)

random_field <- gp(x = coords, # lon, lat, year
                   kernel = kernel,
                   inducing = unique(coords)  # shouldn't do much to simulated data but might be very helpful for actual data
)

gp_mean_obs <- X_obs %*% beta + random_field

mut_prob_obs <- ilogit(gp_mean_obs)

distribution(mut_data$present) <- greta::binomial(size = mut_data$tested,
                                                  prob = mut_prob_obs)

# beta = parameters$beta
# m <- model(parameters$kernel_lengthscale_space, 
#            parameters$kernel_lengthscale_time, 
#            parameters$kernel_sd, 
#            parameters$beta)
m <- model(kernel_lengthscale_space, 
           kernel_lengthscale_time, 
           kernel_sd, 
           beta)

# get a fun little error down here ...
draws <- mcmc(m, n_samples = 3000, 
               initial_values = replicate(4, 
                                          initials(beta = rep(0,3)), 
                                          simplify = FALSE))

calculate(mut_prob_obs, values = list(kernel_lengthscale_space = 1,
                                       kernel_lengthscale_time = 1,
                                       kernel_sd = 0.1,
                                       beta = c(0, 1, 1)))

# alpha <- normal(0, 1)
# beta <- normal(0, 1)
# sigma <- lognormal(1, 0.1)
# y <- as_data(iris$Petal.Width)
# mu <- alpha + iris$Petal.Length * beta
# distribution(y) <- normal(mu, sigma)
# m <- model(alpha, beta, sigma)
# 
# # sample values of the parameters, or different observation data (y), from
# # the priors (useful for prior # predictive checking) - see also
# # ?simulate.greta_model
# calculate(alpha, beta, sigma, nsim = 100)
# calculate(y, nsim = 100)
# 
# # calculate intermediate greta arrays, given some parameter values (useful
# # for debugging models)
# calculate(mu[1:5], values = list(alpha = 1, beta = 2, sigma = 0.5))
# calculate(mu[1:5], values = list(alpha = -1, beta = 0.2, sigma = 0.5))

# Error in `self$check_reasonable_starting_values()`:
# ! Could not find reasonable starting values after 20 attempts.
# Please specify initial values manually via the `initial_values` argument
# suggests everything is so bad it won't be accepted
# get calculate working and test what gets crapped out ..
# `valid_parameters`

png("output/circmat_trace.png", height=750, width=1500)
bayesplot::mcmc_trace(draws)
dev.off()

r_hats <- coda::gelman.diag(draws,
                            autoburnin = FALSE,
                            multivariate = FALSE)
summary(r_hats$psrf)

# save everything and do the prediction in a separate script
write_rds(mut_data, "output/circmat_model/mut_data.rds")
write_rds(parameters, "output/circmat_model/parameters.rds")
# write_rds(mut_prob_obs, "output/circmat_model/mut_prob_obs.rds")
write_rds(kernel, "output/circmat_model/kernel.rds")
write_rds(random_field, "output/circmat_model/random_field.rds")
write_rds(m, "output/circmat_model/m.rds")
write_rds(draws, "output/circmat_model/draws.rds")




# would be nice to have a figure that shows difference between absolute distance 
# and great circle

pts = cbind(rep(20, 10), seq(30, -30, length.out = 10))
dists = data.frame(euclid = as.matrix(dist(pts))[1,],
                   gc = distHaversine(pts, c(20, 30))) %>%
  mutate(euclid = euclid / , # this is in degrees ...
         gc = gc/1000) # km

plot(0:9, dists[,1])
lines(0:9, dists[,2])

# okay so this is a bust so far ...
