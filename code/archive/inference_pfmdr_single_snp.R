# Fit pfmdr86 model and write outputs to output/circmat_pfmdr86/directory
source("code/setup.R")
source("code/build_design_matrix.R")

snp <- 1246

mut_data <- setup_mut_data(paste0("../moldm/clean/pfmdr_single_", snp, ".csv"))

out <- build_design_matrix(covariates,
                           coords = mut_data,
                           scale = FALSE,
                           temporal_var = TRUE,
                           temporal_range = pfpr_years,
                           degs_to_rads = TRUE)
X_obs <- out$df
scaled_years <- out$scaled_years

coord_cols <- c("x_rd", "y_rd", "year_scaled")

design_cols <- c("intercept", "year_scaled", "pfpr")

# hyperparameters
circmat_len <- greta::lognormal(meanlog = -2, sdlog = 1)
circmat_var <- greta::lognormal(meanlog = -2, sdlog = 1)
expo_len <- normal(0, 1, truncation = c(0, Inf))
expo_var <- normal(0, 1, truncation = c(0, Inf))
nugget_var <- normal(0,1, truncation = c(0, Inf))

# kernel & GP
kernel <- circmat(circmat_len, variance = circmat_var, columns = 1:2) + 
  rbf(lengthscales = expo_len, variance = expo_var, columns = 3) +
  constant(nugget_var)

kmn <- kmeans(X_obs[,coord_cols], centers = 40)
random_field <- gp(x = X_obs[,coord_cols], 
                   kernel = kernel,
                   inducing = kmn$centers)

beta <- normal(0, 1, dim = 3)
gp_mean_obs <- X_obs[,design_cols] %*% beta + random_field
X_prob_obs <- ilogit(gp_mean_obs)

# likelihood
distribution(X_obs$present) <- binomial(X_obs$tested, X_prob_obs)

# fit the model by Hamiltonian Monte Carlo
m <- model(circmat_len, circmat_var, expo_len, expo_var, nugget_var, beta)

# if the outcome is this sensitive to the seed I'm using this is cause for alarm ...
# make sure you use 141 for pfmdr184 ....
# I might need to reconfigure model setup here ....
set.seed(141)
# require a different seed to get pfmdr1246 to work
set.seed(1246)
draws <- mcmc(m, 
              n_samples = 3000,
              initial_values = initials(circmat_len = 0.02,
                                        circmat_var = 0.5,
                                        expo_var = 0.5,
                                        expo_len = 1,
                                        nugget_var = 0.5,
                                        beta = rep(0, 3)))
draws <- extra_samples(draws, 3000)
draws <- extra_samples(draws, 3000)

bayesplot::mcmc_trace(draws)

r_hats <- coda::gelman.diag(draws,
                            autoburnin = FALSE,
                            multivariate = FALSE)
summary(r_hats$psrf)

# identifiability an issue between temporal lengthscale and variance

parameters <- list(circmat_len, circmat_var, expo_len, expo_var, nugget_var, beta)
names(parameters) <- c("circmat_len", "circmat_var", "expo_len", "expo_var",
                       "nugget_var", "beta")

# save everything and do the prediction in a separate script
write_rds(mut_data, paste0("output/circmat_pfmdr", snp, "/mut_data.rds"))
write_rds(parameters, paste0("output/circmat_pfmdr", snp, "/parameters.rds"))
write_rds(kernel, paste0("output/circmat_pfmdr", snp, "/kernel.rds"))
write_rds(random_field, paste0("output/circmat_pfmdr", snp, "/random_field.rds"))
write_rds(m, paste0("output/circmat_pfmdr", snp, "/m.rds"))
write_rds(draws, paste0("output/circmat_pfmdr", snp, "/draws.rds"))

# quick little visualisation + prediction script to one year:
# ras_agg <- covariates$pfpr_2019 %>%
#   aggregate(fact = 10)
# # could write this into `build_design_matrix` ...
# X_pixel <- ras_agg %>%
#   xyFromCell(terra::cells(ras_agg)) %>%
#   as.data.frame() %>%
#   mutate(pfpr = unlist(extract(ras_agg, terra::cells(ras_agg))),
#          year = 2019,
#          year_scaled = scaled_years$`2019`,
#          x_rd = degrees_to_radians(x),
#          y_rd = degrees_to_radians(y),
#          intercept = 1)
# ras_agg_rd <- ras_agg
# ext(ras_agg_rd) <- degrees_to_radians(as.vector(ext(ras_agg)))
# 
# plot(ras_agg)
# points(X_obs[X_obs$present == 0,c("x", "y")], col="red", cex=0.5)
# points(X_obs[,c("x", "y")], cex = X_obs$present / X_obs$tested * 20, col="orange")
# 
# plot(ras_agg_rd)
# points(X_obs[X_obs$present == 0,c("x_rd", "y_rd")], col="red", cex=0.5)
# points(X_obs[,c("x_rd", "y_rd")], cex = X_obs$present / X_obs$tested * 20, col="orange")
# f_plot <- greta.gp::project(random_field, X_pixel[,coord_cols])
# mut_freq_plot <- X_pixel[,design_cols] %*% beta + f_plot
# y_plot <- greta::calculate(mut_freq_plot,
#                            values = draws)
# med_vals <- apply(y_plot[[1]], 2, median)
# sd_vals <- apply(y_plot[[1]], 2, sd)
# 
# med_ras <- ras_agg
# med_ras[cells(med_ras)] <- unlist(calculate(ilogit(med_vals)))
# sd_ras <- ras_agg
# sd_ras[cells(sd_ras)] <- sd_vals
# 
# par(mfrow=c(1,2))
# plot(med_ras, main="Median posterior samp")
# points(X_obs[X_obs$present == 0,c("x", "y")], col="red", cex=0.5)
# points(X_obs[,c("x", "y")], cex = X_obs$present / X_obs$tested * 20, col="orange")
# 
# plot(sd_ras, main="SD posterior samp")
# points(X_obs[X_obs$present == 0,c("x", "y")], col="red", cex=0.5)
# points(X_obs[,c("x", "y")], cex = X_obs$present / X_obs$tested * 20, col="orange")


