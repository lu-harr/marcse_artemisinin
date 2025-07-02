# script for multiple binomial regression ...?
setwd("~/Desktop/MARCSE/k13_seafrica")
source("code/setup.R")
source("code/build_design_matrix.R")

# aspirational:
# distribution(snps$positive) <- binomial(snps$tested,)

mut_data <- setup_multiple_snps("../moldm/clean/pfmdr_single_locus.csv")

out <- build_design_matrix(covariates,
                           coords = mut_data,
                           scale = FALSE,
                           temporal_var = TRUE,
                           temporal_range = pfpr_years,
                           degs_to_rads = TRUE)
X_obs <- out$df
scaled_years <- out$scaled_years

n_snp <- length(unique(mut_data$snp))
n_latent <- 2  # guessing ...
coord_cols <- c("x_rd", "y_rd")#, "year_scaled")
design_cols <- c("intercept", "year_scaled", "pfpr")

source("code/define_greta_params.R")

# tmp <- pfpr$pfpr_2000
# values(tmp)[!is.na(values(tmp))] <- cells(tmp)
# mut_data$coord_id <- unlist(terra::extract(tmp, mut_data[, c("x", "y")], ID = FALSE))
# coord_id is in fact the rownumber of the entry in mut_data in coords I think

# where do I encode the number of betas?
parameters <- define_greta_parameters(n_snp = n_snp,
                                      n_latent = n_latent)

# matern 5/2 isotropic kernel - space-only
kernel_lengthscale <- normal(14, 1, truncation = c(0, Inf))
kernel_sd <- normal(0, 1, truncation = c(0, Inf))
kernel <- mat52(lengthscales = c(kernel_lengthscale, kernel_lengthscale),
                variance = kernel_sd ^ 2)

coords <- X_obs[,coord_cols] %>%
  mutate(coord_id = row_number())

# define knots for reduced-rank GP approximation
kmn <- kmeans(coords, centers = 25)

# define GPs over spatial latent factors, evaluated at all data locations
latents_obs <- gp(x = coords,
                  kernel = kernel,
                  inducing = kmn$centers,
                  n = n_latent)

# extract out the design matrices (pre-scaled)
# X_obs <- build_design_matrix(X_obs[,design_cols], X_obs[,coord_cols])

# combine these with parameters to get matrices SNP frequencies at the SNP data locations
snp_freq_logit_obs <- X_obs[,design_cols] %*% t(parameters$beta) +
  t(parameters$loadings %*% t(latents_obs))
snp_freq_obs <- ilogit(snp_freq_logit_obs)
# what is the dimension of this thing? should be nrow(X_obs) * ncol(beta)

# define the likelihood over the SNPs
snp_data_index <- cbind(X_obs$coord_id, mut_data$snp_id)
distribution(mut_data$present) <- greta::binomial(size = mut_data$tested,
                                                   prob = snp_freq_obs[snp_data_index]) # is there a way to unpack this?

beta <- parameters$beta
m <- model(kernel_lengthscale, kernel_sd, beta)
draws <- mcmc(m)
draws <- extra_samples(draws, 2000)

# to my great surprise this appears to run
r_hats <- coda::gelman.diag(draws,
                            autoburnin = FALSE,
                            multivariate = FALSE)
summary(r_hats$psrf)











