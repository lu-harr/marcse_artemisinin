#setwd("~/Desktop/MARCSE/k13_seafrica")
source("code/setup.R")

coarse_mat <- aggregate(pfpr$pfpr_2000, fact = 10)

# or install local version ...:
#devtools::install_github("lu-harr/greta.gp.st.on.earth")
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
#     # JF/YSF use flat priors - could make these flatter?
#     beta = normal(0, 3, dim=3), # I have three betas so dim is 3 (intercept, time, PfPR)
#     lengthscale_space = lognormal(-0.5, 1),
#     # YSF's prior in degrees:
#     # degrees_to_radians(qlnorm(c(0.005, 0.995), meanlog = 0.5, sdlog = 1)) # [0.13, 21.67] degs
#     # This prior in radians is in a similar order of magnitude:
#     # qlnorm(c(0.005, 0.995), meanlog = -3.5, sdlog = 1)
#     
#     lengthscale_time = lognormal(-0.5, 1),
#     # [0.5%, 99.5%] quantile of [0.046, 7.971] years
#     # qlnorm(c(0.005, 0.995), meanlog = -0.5, sdlog = 1)
#     
#     kernel_sd = normal(0, 2, truncation = c(0, Inf))
#     # YSF uses HalfNormal (Normal truncated at 0)
#   )
# }

# parameters <- define_greta_parameters() # are there any choices that should be arguments?

# JF/YSF use flat priors - could make these flatter?
beta = normal(0, 3, dim=3) # I have three betas so dim is 3 (intercept, time, PfPR)
lengthscale_space = lognormal(meanlog = -3, sdlog = 1.5)
# YSF's prior in degrees:
# degrees_to_radians(qlnorm(c(0.005, 0.995), meanlog = 0.5, sdlog = 1)) # [0.13, 21.67] degs
# This prior in radians is in a similar order of magnitude:
# qlnorm(c(0.005, 0.995), meanlog = -3.5, sdlog = 1)
qlnorm(c(0.005, 0.995), meanlog = -0.5, sdlog = 1.5)

lengthscale_time = lognormal(-3, 1)
# [0.5%, 99.5%] quantile of [0.046, 7.971] years
# temporal correlation is 0.05 when delta = 6*lengthscale
# hist(X_obs[,2], breaks=20)
# qlnorm(c(0.005, 0.995), meanlog = -0.5, sdlog = 1) # (unscaled years)
# abline(v = median(X_obs[,2]), lwd = 2)
# abline(v = median(X_obs[,2]) + c(1,-1)*qlnorm(0.005, meanlog = -3, sdlog = 1)*6, col="red")
# abline(v = median(X_obs[,2]) + c(1,-1)*qlnorm(0.9, meanlog = -3, sdlog = 1)*6, col="blue")
# qlnorm(c(0.005, 0.995), meanlog = -3, sdlog = 1)

kernel_sd = normal(0, 1, truncation = c(0, Inf))
# YSF uses HalfNormal (Normal truncated at 0)

# kernel <- circmat(lengthscale = parameters$lengthscale_space,
#                   variance = 1,
#                   columns = c(1, 2),
#                   circumference = 1) +
#           expo(lengthscales = parameters$lengthscale_time,
#                variance = 1,
#                columns = 3) +
#   white(parameters$kernel_sd**2)

kernel <- circmat(lengthscale = lengthscale_space,
                  variance = 1,
                  columns = c(1, 2),
                  circumference = 1) +
  expo(lengthscales = lengthscale_time,
       variance = 1,
       columns = 3) +
  white(kernel_sd**2)

# need to incorporate `kernel _sd` .. could go in both variance spots ..?

coords <- coords %>%
  dplyr::select(x, y, year)

#kmn <- kmeans(coords, centers = 15)

random_field <- gp(x = coords, # lon, lat, year
                   kernel = kernel,
                   inducing = unique(coords))#,  # shouldn't do much to simulated data but might be very helpful for actual data
                   #inducing = kmn$centers)

# gp_mean_obs <- X_obs %*% parameters$beta + random_field
gp_mean_obs <- X_obs %*% beta + random_field

mut_prob_obs <- ilogit(gp_mean_obs)

distribution(mut_data$present) <- greta::binomial(size = mut_data$tested,
                                                  prob = mut_prob_obs)

# m <- model(parameters$lengthscale_space,
#            parameters$lengthscale_time,
#            parameters$kernel_sd,
#            parameters$beta)
m <- model(lengthscale_space,
           lengthscale_time,
           kernel_sd,
           beta)

# draws <- mcmc(m, n_samples = 100, warmup = 100, 
#               initial_values = replicate(4, initials(`parameters$beta` = rep(0,3),
#                                                      `parameters$lengthscale_space` = 0.01,
#                                                      `parameters$lengthscale_time` = 0.01,
#                                                      `parameters$kernel_sd` = 0.5), simplify = FALSE))
draws <- mcmc(m, n_samples = 100, warmup = 100, 
              initial_values = replicate(4, initials(`beta` = rep(0,3),
                                                     `lengthscale_space` = 1,
                                                     `lengthscale_time` = 1,
                                                     `kernel_sd` = 0.1), simplify = FALSE))

coords_pixel <- coarse_mat %>%
  terra::xyFromCell(cell = terra::cells(coarse_mat)) %>%
  as.data.frame() %>%
  degrees_to_radians() %>%
  mutate(scaled_year = 1.0320937)

f_plot <- project(random_field, coords_pixel)
# Make sure you populate nsim here:

prior_eg = greta::calculate(f_plot, values = list(lengthscale_space = 2,
                                                        lengthscale_time = 2,
                                                        kernel_sd = 2,
                                                        beta = c(0, 0, 0)),
                            nsim = 1)
tmp = coarse_mat
tmp[!is.na(coarse_mat)] <- unlist(prior_eg)

# with much lower lengthscale:
plot(tmp,
     xlab = "Longitude", ylab = "Latitude")
# the new palette is sending me
#points(x_rd, cex = z - min(z) + 1)

#png("output/circmat_trace.png", height=750, width=1500)
bayesplot::mcmc_trace(draws)
#dev.off()

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
