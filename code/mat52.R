source("code/build_design_matrix.R")
source("code/setup.R")

message("test")

#install_greta_deps()

#message("test1")

library(greta.gp) # we won't need my diy version for this one

# a start - works okay
coords <- mut_data %>%
  dplyr::select(x, y, year, scaled_year)

source("code/build_design_matrix.R")
X_obs <- build_design_matrix(covariates, 
                             coords, 
                             scale = FALSE, 
                             temporal_range = pfpr_years)
# we want this to give us back scaled years but need to provide unscaled years
# so that it can grab from correct annual covt raster
message(nrow(X_obs))


define_greta_parameters <- function(){
  list(
    beta = normal(0, 3, dim=3), # I have three betas so dim is 3 (intercept, time, PfPR)
    kernel_lengthscale_space = normal(0, 3, truncation = c(0, Inf)),
    kernel_lengthscale_time = normal(0, 3, truncation = c(0, Inf)),
    kernel_sd = normal(0, 1, truncation = c(0, Inf))
  )
}

parameters <- define_greta_parameters() # are there any choices that should be arguments?
# beta = normal(0, 3, dim = 3)
# kernel_lengthscale_space = normal(0, 3, truncation = c(0, Inf))
# kernel_lengthscale_time = normal(0, 3, truncation = c(0, Inf))
# kernel_sd = normal(0, 1, truncation = c(0, Inf))


kernel <- mat52(lengthscales = c(parameters$kernel_lengthscale_space, 
                                 parameters$kernel_lengthscale_space, 
                                 parameters$kernel_lengthscale_time),
                variance = parameters$kernel_sd ^ 2)

coords <- coords %>%
  dplyr::select(x, y, scaled_year)

# dhps data seems to need inducing points:
#kmn <- kmeans(coords, centers = 15)

random_field <- gp(x = coords, # lon, lat, year
                   kernel = kernel,
                   #inducing = kmn$centers
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


library(GGally)
library(brms)

modified_density = function(data, mapping, ...) {
  ggally_densityDiag(data, mapping, ...) + 
    scale_fill_brewer(type = "qual", palette = "Purples")
}

# ggplot is the silliest darn software going
# why would anyone want or need to change a colour palette?
# beats me
modified_points = function(data, mapping, ...) {
  ggally_points(data, mapping, ...) + 
    scale_color_brewer(palette = "Purples")
}

modified_cor = function(data, mapping, ...){
  ggally_cor(data, mapping, ...) +
    scale_color_brewer(palette = "Purples")
}

post <- as_draws_df(draws) %>%
  rename("Lengthscale (spatial)" = "parameters$kernel_lengthscale_space",
         "Lengthscale (temporal)" = "parameters$kernel_lengthscale_time",
         "sigma" = "parameters$kernel_sd",
         "beta_0" = "parameters$beta[1,1]",
         "beta_1" = "parameters$beta[2,1]",
         "beta_2" = "parameters$beta[3,1]")

library(bayesplot)
color_scheme_set("brewer-Purples")
bayesplot::mcmc_trace(post)
ggsave("figures/mat52_trace.png", height=10, width=15)

# love the purps but they're a bit too light to use for cor
post %>%
  mutate(chain = as.factor(.chain)) %>%
  ggpairs(columns = 1:6,
          mapping = aes(colour = chain),
          lower = list(continuous = wrap(modified_points, alpha = 0.4)), # can't seem to turn alpha down here?
          diag = list(continuous = wrap(modified_density, alpha = 0.5)),
          upper = list(continuous = modified_cor)) # probably don't need corrs for chains?
ggsave("figures/chain_corr_mat52.png", height=12, width=12)

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

# I've seen worse .. just needs to run a bit longer


# # save everything and do the prediction in a separate script
# write_rds(mut_data, "../dhps_mapping/quintuple_question/quintuple/mut_data.rds")
# write_rds(parameters, "../dhps_mapping/quintuple_question/quintuple/parameters.rds")
# # write_rds(mut_prob_obs, "../dhps_mapping/quintuple_question/quintuple/mut_prob_obs.rds")
# write_rds(kernel, "../dhps_mapping/quintuple_question/quintuple/kernel.rds")
# write_rds(random_field, "../dhps_mapping/quintuple_question/quintuple/random_field.rds")
# write_rds(m, "../dhps_mapping/quintuple_question/quintuple/m.rds")
# write_rds(draws, "../dhps_mapping/quintuple_question/quintuple/draws.rds")
# 
# # save everything and do the prediction in a separate script
# write_rds(mut_data, "../dhps_mapping/quintuple_question/540E/mut_data.rds")
# write_rds(parameters, "../dhps_mapping/quintuple_question/540E/parameters.rds")
# # write_rds(mut_prob_obs, "../dhps_mapping/quintuple_question/540E/mut_prob_obs.rds")
# write_rds(kernel, "../dhps_mapping/quintuple_question/540E/kernel.rds")
# write_rds(random_field, "../dhps_mapping/quintuple_question/540E/random_field.rds")
# write_rds(m, "../dhps_mapping/quintuple_question/540E/m.rds")
# write_rds(draws, "../dhps_mapping/quintuple_question/540E/draws.rds")






