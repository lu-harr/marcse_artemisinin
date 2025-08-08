# gneiting inference slurm - aiming for same model, so should be able to provide
# marker to command line
start <- Sys.time()

source("code/setup.R")
source("code/build_design_matrix.R")
source("code/predict_to_raster.R")

snp = "k13"
seed = 123

out_dir <- paste0(snp, "/gneiting_ahmc/")

in_dat <- ifelse(snp == "k13",
                 "data/clean/moldm_k13_nomarker.csv",
                 ifelse(snp == "crt76",
                        "data/clean/moldm_crt76.csv",
                        ifelse(snp == "k13_marcse",
                               "data/clean/moldm_marcse_k13_nomarker.csv",
                               paste0(paste0("data/clean/pfmdr_single_", snp, ".csv")))))


message(paste0("Reading in from: ", in_dat))
message(getwd())
message("Enforcing min year for surveyor data - 2000")
mut_data <- setup_mut_data(in_dat, min_year = 2000)
mut_data <- mut_data %>% filter(y <=0)

out <- build_design_matrix(covariates,
                           coords = mut_data,
                           scale = FALSE,
                           temporal_var = TRUE,
                           temporal_covt_range = pfpr_years,
                           degs_to_rads = TRUE)
X_obs <- out$df
scaled_years <- out$scaled_years
message(scaled_years)

coord_cols <- c("x_rd", "y_rd", "year_scaled")

design_cols <- c("intercept", "year_scaled", "pfpr")

# hyperparameters
gneiting_len <- greta::normal(0, 3, truncation = c(0, Inf))
gneiting_tim <- greta::normal(0, 3, truncation = c(0, Inf))
gneiting_sd <- greta::normal(0, 2, truncation = c(0, Inf))
nugget_sd <- greta::normal(0, 3, truncation = c(0, Inf)) # Median :1.041  Mean   :1.115

# gneiting_len <- greta::lognormal(0, 1)
# gneiting_tim <- greta::lognormal(0, 1)
# gneiting_sd <- greta::lognormal(0, 1)
# nugget_sd <- greta::lognormal(0, 1)

# kernel & GP
# could potentially give stricter priors to variance/nugget here - trouble with IDability?
kernel <- greta.gp::gneiting(lengthscale = gneiting_len, 
                   timescale = gneiting_tim,
                   variance = gneiting_sd ** 2, 
                   columns = 1:3) + 
  white(nugget_sd ** 2)

kmn <- kmeans(X_obs[,coord_cols], centers = 40)
random_field <- greta.gp::gp(x = X_obs[,coord_cols], 
                   kernel = kernel,
                   inducing = kmn$centers)

beta <- greta::normal(0, 1, dim = 3)
gp_mean_obs <- X_obs[,design_cols] %*% beta + random_field
X_prob_obs <- ilogit(gp_mean_obs)

# likelihood
distribution(X_obs$present) <- binomial(X_obs$tested, X_prob_obs)

# fit the model by Hamiltonian Monte Carlo
m <- model(gneiting_len, gneiting_tim, gneiting_sd, nugget_sd, beta)

set.seed(seed)
draws <- mcmc(m,
              n_samples = 3000,
              initial_values = initials(gneiting_len = 1,
                                        gneiting_tim = 3,
                                        gneiting_sd = 13,
                                        nugget_sd = 0.5,
                                        beta = rep(0, 3)))


#draws <- extra_samples(draws, 20000)

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

end <- Sys.time()

print(end - start)


AGG_FACTOR = 1

stable_transmission_mask <- rast("data/stable_transmission_mask.grd") %>%
  aggregate(AGG_FACTOR)

for (year in c(2006, 2010, 2014, 2018, 2022)){
  preds <- predict_to_ras(covariates,
                          year,
                          draws,
                          parameters,
                          random_field,
                          agg_factor = AGG_FACTOR,
                          scaled_year = scaled_years[[as.character(year)]],
                          coord_cols = c("x_rd", "y_rd", "year_scaled"),
                          design_cols = c("intercept", "year_scaled", "pfpr"),
                          stable_transmission_mask = stable_transmission_mask)
  
  writeRaster(preds, paste0("output/", out_dir, year, "_preds.grd"), overwrite = TRUE)
}




#message(out_dir)
out_dir <- paste0("output/k13/gneiting_ahmc/")
to_read <- grep("^20.*\\.grd$", list.files(out_dir), value = TRUE)
#message(to_read)
writeRaster(x = rast(paste0(out_dir, to_read)), 
            filename = paste0(out_dir, "preds_all.grd"), overwrite = TRUE)



preds <- rast("output/k13/gneiting_ahmc/preds_all.grd")
# yeah I mean says it all really
plot(preds$`2006_post_median`)
plot(preds$`2006_post_sd`)
plot(preds$`2010_post_median`)
plot(preds$`2010_post_sd`)
plot(preds$`2014_post_median`)
plot(preds$`2014_post_sd`)
plot(preds$`2018_post_median`)
plot(preds$`2018_post_sd`)
plot(preds$`2022_post_median`)
plot(preds$`2022_post_sd`)

# play around with betabinomial??
# what does var look like if I turn up rho?
source("code/betabinomial_p_rho.R")
# another parameterisation in brms ??

p <- seq(0, 100)
n <- 100

plot(p, dbbinom(p, n, 10, 1), type = "l")
lines(p, dbbinom(p, n, 1, 10))
lines(p, dbbinom(p, n, 0.5, 0.5))
lines(p, dbbinom(p, n, 2,2))
lines(p, dbbinom(p, n, 10, 10))










