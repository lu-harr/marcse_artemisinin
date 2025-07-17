# Fit K13 model and write outputs to output/gneiting_k13/directory
#setwd("~/Desktop/MARCSE/k13_seafrica")
source("code/setup.R")
source("code/build_design_matrix.R")

out_dir <- "k13/gneiting_sparse/"

mut_data <- setup_mut_data("data/clean/moldm_k13_nomarker.csv", min_year = 2000)

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
gneiting_len <- normal(0, 3, truncation = c(0, Inf))
gneiting_tim <- normal(0, 3, truncation = c(0, Inf))
gneiting_sd <- normal(0, 2, truncation = c(0, Inf))
nugget_sd <- normal(0, 3, truncation = c(0, Inf))

# kernel & GP
kernel <- gneiting(lengthscale = gneiting_len, 
                   timescale = gneiting_tim,
                   variance = gneiting_sd ** 2, 
                   columns = 1:3) + 
  white(nugget_sd ** 2)

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
m <- model(gneiting_len, gneiting_tim, gneiting_sd, nugget_sd, beta)

set.seed(123)
draws <- mcmc(m,
              n_samples = 10000,
              initial_values = initials(gneiting_len = 1,
                                        gneiting_tim = 3,
                                        gneiting_sd = 13,
                                        nugget_sd = 0.5,
                                        beta = c(-0.8, 1.5, 0.5)))


draws <- extra_samples(draws, 10000)

#bayesplot::mcmc_trace(draws)

r_hats <- coda::gelman.diag(draws,
                            autoburnin = FALSE,
                            multivariate = FALSE)
summary(r_hats$psrf)

parameters <- list(gneiting_len, gneiting_tim, gneiting_sd, nugget_sd, beta)
names(parameters) <- c("gneiting_len", "gneiting_tim", "gneiting_sd",
                       "nugget_sd", "beta")

# save everything and do the prediction in a separate script
write_rds(mut_data, paste0("output/", out_dir, "mut_data.rds"))
write_rds(parameters, paste0("output/", out_dir, "parameters.rds"))
write_rds(kernel, paste0("output/", out_dir, "kernel.rds"))
write_rds(random_field, paste0("output/", out_dir, "random_field.rds"))
write_rds(m, paste0("output/", out_dir, "m.rds"))
write_rds(draws, paste0("output/", out_dir, "draws.rds"))

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


