#' Make model predictions to a raster of pixels
#'
#' @param stack Rast of covariates, with names containing years
#' @param year numeric
#' @param draws greta_mcmc_list of simulations from model posterior
#' @param parameters list of greta arrays, defined in inference step
#' @param random_field greta array, defined in inference step
#' @param scaled_year numeric: `year` with scaling applied
#' @param agg_factor numeric: aggregate `stack` before predicting to it? If no, set to 1.
#' @param nsim numeric: number of simulations to draw from predictions
#' @param stable_transmission_mask Rast or `NULL`: If Rast provided, predictions 
#' are masked to areas of stable malaria transmission.
#' @param coord_cols vector of strings: names of coordinate columns in design matrix
#' @param design_cols vector of strings: names of covariate columns in design matrix
#'
#' @returns Rast
#' @export
#'
#' @examples
#' \dontrun{
#' predict_to_ras(covariates, 
#'                2010, 
#'                draws, 
#'                parameters,
#'                random_field,
#'                agg_factor = 15,
#'                stable_transmission_mask = stable_transmission_mask)
#' }
#' 
predict_to_ras <- function(stack, 
                           year, 
                           draws, 
                           parameters,
                           random_field,
                           scaled_year = 0,
                           agg_factor = 1,
                           nsim = 100,
                           stable_transmission_mask = NULL,
                           coord_cols = c("x", "y", "year"),
                           design_cols = c("year")){
  cov_years <- names(stack)
  cov_years <- as.numeric(gsub("[^0-9]", "", cov_years))
  lab_year <- year
  if (year > max(cov_years)){
    year = max(cov_years)
  } else if (year < min(cov_years)){
    year = min(cov_years)
  }
  # message(year)
  
  # retrieve raster for `year`
  ras <- stack[[grep(as.character(year), names(stack))]]
  
  # aggregate covariate raster for `year`
  if(agg_factor != 1){ras <- aggregate(ras, fact = agg_factor)}
  
  # make up df of coordinates and covariate values
  coords <- cbind(terra::xyFromCell(ras, cell = terra::cells(ras)),
                  rep(year, length(terra::cells(ras)))) %>%
    as.data.frame()
  
  names(coords) <- c("x", "y", "year") # this is *not* coord_cols
  # message(names(coords))
  
  tmp <- build_design_matrix(ras,
                             coords,
                             temporal_var = FALSE,
                             scale = FALSE,
                             degs_to_rads = TRUE)
  
  # finish off design matrix
  X_pixel <- tmp$df %>%
    dplyr::mutate(year_scaled = scaled_year)
  message(paste0("Scaled year: ", scaled_year))
  # X_pixel <- dplyr::mutate(X_pixel, year_scaled = scaled_year)
  message(coord_cols) # this is xyyear
  message(design_cols) # this is interceptyear_scaledpfpr
  message(paste(names(X_pixel), collapse = ", "))
  message(dim(X_pixel[,coord_cols]))
  message(dim(X_pixel[,design_cols]))
  # 
  # project random field to coordinates we would like predictions for
  random_field_pixel <- greta.gp::project(random_field, X_pixel[,coord_cols])
  
  # message(dim(t(parameters$beta)))
  # message(dim(random_field_pixel))
  mut_freq_pixel <- (X_pixel[,design_cols] %*% parameters$beta + 
                       # transform here? Is R taking over?
                       random_field_pixel) %>%
    ilogit()
  # message(dim(mut_freq_pixel))
  
  post_pixel_sims <- greta::calculate(mut_freq_pixel,
                                      values = draws,
                                      nsim = 500,
                                      trace_batch_size = 1) # reducing: will take longer, use less mem
  # message("now here")
  probs = c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  post_pixel_quants = apply(post_pixel_sims$mut_freq_pixel[,,1], 2, 
                             quantile, 
                             probs = probs)
  post_pixel_mean <- apply(post_pixel_sims$mut_freq_pixel[, , 1], 2, mean)
  post_pixel_sd <- apply(post_pixel_sims$mut_freq_pixel[, , 1], 2, sd)
  post_pixel_sd_logit <- apply(post_pixel_sims$mut_freq_pixel[,,1], 2, function(x){sd(log(x / (1 - x)))})
  post_summary <- cbind(t(post_pixel_quants), post_pixel_mean, post_pixel_sd, post_pixel_sd_logit)
  
  out <- ras * rep(0, nrow(post_pixel_quants) + 3) # assuming we have at least two layers in there ..
  names(out) <- paste0(lab_year, "_", c(probs*100, "mean", "sd", "sdscaled"))
  out[terra::cells(out)] <- post_summary
  
  if(!is.null(stable_transmission_mask)){
  # if (length(unique(suppressWarnings(values(stable_transmission_mask)))) != 1){
    # let it be known that I did some googling about this :(
    # in raster I would have plopped `raster(NA)` in the function definition
    out <- mask(out, stable_transmission_mask)
  }
  out
}

# out_dir <- "output/crt76/gneiting_sparse/"
# source("code/setup.R")
# source("code/build_design_matrix.R")
# scaled_years <- scale_years(range(pfpr_years))
# 
# # bring in all of the other outputs here too
# AGG_FACTOR <- 5
# mut_data <- read_rds(paste0(out_dir, "mut_data.rds"))
# stable_transmission_mask <- rast("data/stable_transmission_mask.grd") %>%
#   aggregate(AGG_FACTOR)
# random_field <- read_rds(paste0(out_dir, "random_field.rds"))
# parameters <- read_rds(paste0(out_dir, "parameters.rds"))
# draws <- read_rds(paste0(out_dir, "draws.rds"))
# 
# tmp <- predict_to_ras(covariates,
#                        2023,
#                        draws,
#                        parameters,
#                        random_field,
#                        agg_factor = AGG_FACTOR,
#                        stable_transmission_mask = stable_transmission_mask,
#                       design_cols = c("intercept", "year_scaled", "pfpr"))
# 
# plot(tmp)
# 
# library(tidyterra)
# ggplot() +
#   geom_spatraster(data = tmp) +
#   scale_fill_distiller(palette = "RdBu") +
#   facet_wrap(~lyr)
# 
# library(looseVis)
# looseVis::rast_plot(tmp$`2023_0`)


predict_to_ras_hier <- function(stack, 
                           year, 
                           draws, 
                           parameters,
                           random_field,
                           scaled_year = 0,
                           agg_factor = 1,
                           nsim = 100,
                           stable_transmission_mask = NULL,
                           coord_cols = c("x", "y", "year"),
                           design_cols = c("year"),
                           n_snp = 1){
  
  ras <- stack[[grep(as.character(year), names(stack))]]
  
  if(agg_factor != 1){ras <- aggregate(ras, fact = agg_factor)}
  
  coords <- cbind(terra::xyFromCell(ras, cell = terra::cells(ras)),
                  rep(year, length(terra::cells(ras)))) %>%
    as.data.frame()
  message(names(coords))
  
  names(coords) <- c("x", "y", "year") # this is *not* coord_cols
  
  tmp <- build_design_matrix(ras,
                             coords,
                             temporal_var = FALSE,
                             scale = FALSE,
                             degs_to_rads = TRUE)
  
  X_pixel <- tmp$df %>%
    mutate(year_scaled = scaled_year)
  
  message(names(X_pixel))
  
  random_field_pixel <- greta.gp::project(random_field, X_pixel[,coord_cols])
  
  # latent factors key distinguishing feature of prediction step:
  mut_freq_pixel <- (X_pixel[,design_cols] %*% t(parameters$beta) + 
                     t(parameters$loadings %*% t(random_field_pixel))) %>%
    ilogit()
  
  post_pixel_sims <- greta::calculate(mut_freq_pixel,
                                      values = draws,
                                      nsim = 100,
                                      trace_batch_size = 50) # reducing: will take longer, use less mem
  
  out <- rep(ras, n_snp * 2) * 0
  names(out) <- levels(interaction(year, c("median", "sd"), 1:n_snp))
                       
  for (i in 1:n_snp){ # would be good if I could name snps ...
    post_pixel_median <- apply(post_pixel_sims$mut_freq_pixel[ , , i], 2, median)
    post_pixel_sd <- apply(post_pixel_sims$mut_freq_pixel[ , , i], 2, sd)
    
    out[[i * 2 - 1]][terra::cells(out[[i * 2 - 1]])] <- post_pixel_median
    out[[i * 2]][terra::cells(out[[i * 2]])] <- post_pixel_sd
  }
  
  if(!is.null(stable_transmission_mask)){
    #if (length(unique(suppressWarnings(values(stable_transmission_mask)))) != 1){
    # let it be known that I did some googling about this :(
    # in raster I would have plopped `raster(NA)` in the function definition
    out <- mask(out, stable_transmission_mask)
  }
  out
}


# # for prediction raster:
# AGG_FACTOR = 5
# out_dir <- "output/pfmdr_hier/"
# scaled_years <- scale_years(range(pfpr_years))
# 
# # bring in all of the other outputs here too
# mut_data <- read_rds(paste0(out_dir, "mut_data.rds"))
# stable_transmission_mask <- rast("data/stable_transmission_mask.grd") %>%
#   aggregate(AGG_FACTOR)
# random_field <- read_rds(paste0(out_dir, "latents_obs.rds"))
# parameters <- read_rds(paste0(out_dir, "parameters.rds"))
# draws <- read_rds(paste0(out_dir, "draws.rds"))
# 
# pfpr_years <- c(2010,2015,2020)
# preds <- lapply(pfpr_years, function(year){
#   predict_to_ras_hier(covariates, 
#                       year, 
#                       draws, 
#                       parameters,
#                       random_field,
#                       n_snp = 3,
#                       scaled_year = 0,
#                       agg_factor = AGG_FACTOR,
#                       nsim = 100,
#                       stable_transmission_mask = stable_transmission_mask,
#                       coord_cols = c("x", "y", "year"),
#                       design_cols = c("intercept", "pfpr"))
# })
# 
# preds <- rast(preds)
# plot(preds[[grepl("median", names(preds))]])
# # this all looks like PfPR ....