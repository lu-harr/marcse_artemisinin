source("code/setup.R")
source("code/build_design_matrix.R")
library("greta.gp")

mut_data 

ras_agg <- covariates$pfpr_2019 %>%
  aggregate(fact = 10)

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

# could write this into `build_design_matrix` ...
X_pixel <- ras_agg %>%
  xyFromCell(terra::cells(ras_agg)) %>%
  as.data.frame() %>%
  mutate(pfpr = unlist(extract(ras_agg, terra::cells(ras_agg))),
         year = 2019,
         year_scaled = scaled_years$`2019`,
         x_rd = degrees_to_radians(x),
         y_rd = degrees_to_radians(y),
         intercept = 1)
  
ras_agg_rd <- ras_agg
ext(ras_agg_rd) <- degrees_to_radians(as.vector(ext(ras_agg)))

plot(ras_agg)
points(X_obs[X_obs$present == 0,c("x", "y")], col="red", cex=0.5)
points(X_obs[,c("x", "y")], cex = X_obs$present / X_obs$tested * 20, col="orange")

plot(ras_agg_rd)
points(X_obs[X_obs$present == 0,c("x_rd", "y_rd")], col="red", cex=0.5)
points(X_obs[,c("x_rd", "y_rd")], cex = X_obs$present / X_obs$tested * 20, col="orange")


# hyperparameters
circmat_len <- lognormal(meanlog = -2, sdlog = 1)
circmat_var <- lognormal(meanlog = -2, sdlog = 1)
expo_len <- normal(0, 1, truncation = c(0, Inf))
expo_var <- normal(0, 1, truncation = c(0, Inf))
nugget_var <- normal(0,1, truncation = c(0, Inf))
#obs_sd <- lognormal(0, 2)

# kernel & GP
kernel <- circmat(circmat_len, variance = circmat_var, columns = 1:2) + 
  rbf(lengthscales = expo_len, variance = expo_var, columns = 3) +
  white(nugget_var)

kmn <- kmeans(X_obs[,coord_cols], centers = 15)
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

set.seed(0748)
draws <- mcmc(m, n_samples = 250,
              initial_values = initials(circmat_len = 0.02,
                                        circmat_var = 0.5,
                                        expo_var = 0.5,
                                        expo_len = 0.5,
                                        nugget_var = 0.5,
                                        beta = rep(0, 3)))

bayesplot::mcmc_trace(draws)
# I've seen worse

f_plot <- greta.gp::project(random_field, X_pixel[,coord_cols])
mut_freq_plot <- X_pixel[,design_cols] %*% beta + f_plot
y_plot <- greta::calculate(mut_freq_plot,
                           values = draws)
med_vals <- apply(y_plot[[1]], 2, median)
sd_vals <- apply(y_plot[[1]], 2, sd)

med_ras <- ras_agg
med_ras[cells(med_ras)] <- unlist(calculate(ilogit(med_vals)))
sd_ras <- ras_agg
sd_ras[cells(sd_ras)] <- sd_vals

par(mfrow=c(1,2))
plot(med_ras, main="Median posterior samp")
points(X_obs[X_obs$present == 0,c("x", "y")], col="red", cex=0.5)
points(X_obs[,c("x", "y")], cex = X_obs$present / X_obs$tested * 20, col="orange")

plot(sd_ras, main="SD posterior samp")
points(X_obs[X_obs$present == 0,c("x", "y")], col="red", cex=0.5)
points(X_obs[,c("x", "y")], cex = X_obs$present / X_obs$tested * 20, col="orange")



