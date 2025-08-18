# gneiting inference slurm - aiming for same model, so should be able to provide
# marker to command line
start <- Sys.time()

source("code/setup.R")
source("code/build_design_matrix.R")
source("code/wrap_fit.R")

# this feels a bit unflashy but I can't keep having separate scripts
args <- commandArgs(trailingOnly = TRUE)
snp <- args[1]  # of "k13", "mdr86", "mdr184", "mdr1246", "crt76"
seed <- as.numeric(args[2])
message(paste0("Marker: ", snp))
message(paste0("Seed: ", seed))

# snp = "1246"
# seed = 125

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

fit_binom(mut_data = mut_data,
          covariates = covariates,
          pfpr_years = pfpr_years,
          out_dir = out_dir,
          nchains = 6,
          warmup = 3000,
          nsamples = 20000)

# out <- build_design_matrix(covariates,
#                            coords = mut_data,
#                            scale = FALSE,
#                            temporal_var = TRUE,
#                            temporal_covt_range = pfpr_years,
#                            degs_to_rads = TRUE)
# X_obs <- out$df
# scaled_years <- out$scaled_years
# message(scaled_years)
# 
# coord_cols <- c("x_rd", "y_rd", "year_scaled")
# 
# design_cols <- c("intercept", "year_scaled", "pfpr")
# 
# # hyperparameters
# gneiting_len <- normal(0, 3, truncation = c(0, Inf))
# gneiting_tim <- normal(0, 3, truncation = c(0, Inf))
# gneiting_sd <- normal(0, 2, truncation = c(0, Inf))
# nugget_sd <- normal(0, 3, truncation = c(0, Inf)) # Median :1.041  Mean   :1.115
# 
# # gneiting_len <- greta::lognormal(0, 1)
# # gneiting_tim <- greta::lognormal(0, 1)
# # gneiting_sd <- greta::lognormal(0, 1)
# # nugget_sd <- greta::lognormal(0, 1)
# 
# # kernel & GP
# # could potentially give stricter priors to variance/nugget here - trouble with IDability?
# kernel <- gneiting(lengthscale = gneiting_len, 
#                    timescale = gneiting_tim,
#                    variance = gneiting_sd ** 2, 
#                    columns = 1:3) + 
#   white(nugget_sd ** 2)
# 
# kmn <- kmeans(X_obs[,coord_cols], centers = 40)
# random_field <- greta.gp::gp(x = X_obs[,coord_cols], 
#                    kernel = kernel,
#                    inducing = kmn$centers)
# 
# beta <- normal(0, 1, dim = 3)
# gp_mean_obs <- X_obs[,design_cols] %*% beta + random_field
# X_prob_obs <- ilogit(gp_mean_obs)
# 
# # likelihood
# distribution(X_obs$present) <- binomial(X_obs$tested, X_prob_obs)
# 
# # fit the model by Hamiltonian Monte Carlo
# m <- model(gneiting_len, gneiting_tim, gneiting_sd, nugget_sd, beta)
# 
# {start <- Sys.time()
# set.seed(seed)
# draws <- mcmc(m,
#               sampler = hmc(Lmin = 10, Lmax = 15),
#               chains = 6,
#               warmup = 3000,
#               n_samples = 20000,
#               initial_values = initials(gneiting_len = 1,
#                                         gneiting_tim = 3,
#                                         gneiting_sd = 13,
#                                         nugget_sd = 0.5,
#                                         beta = rep(0, 3)))
# end <- Sys.time()
# end - start}
# 
# 
# draws <- extra_samples(draws, 10000)
# 
# r_hats <- coda::gelman.diag(draws,
#                             autoburnin = FALSE,
#                             multivariate = FALSE)
# summary(r_hats$psrf)
# 
# parameters <- list(gneiting_len, gneiting_tim, gneiting_sd, nugget_sd, beta)
# names(parameters) <- c("gneiting_len", "gneiting_tim", "gneiting_sd",
#                        "nugget_sd", "beta")
# 
# # save everything and do the prediction in a separate script
# write_rds(mut_data, paste0("output/", out_dir, "mut_data.rds"))
# write_rds(parameters, paste0("output/", out_dir, "parameters.rds"))
# write_rds(kernel, paste0("output/", out_dir, "kernel.rds"))
# write_rds(random_field, paste0("output/", out_dir, "random_field.rds"))
# write_rds(m, paste0("output/", out_dir, "m.rds"))
# write_rds(draws, paste0("output/", out_dir, "draws.rds"))

end <- Sys.time()

print(end - start)
message(end - start)