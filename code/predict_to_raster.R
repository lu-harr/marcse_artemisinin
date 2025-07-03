# why no roxygen?

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
  
  ras <- stack[[grep(as.character(year), names(stack))]]
  
  if(agg_factor != 1){ras <- aggregate(ras, fact = agg_factor)}
  
  coords <- cbind(terra::xyFromCell(ras, cell = terra::cells(ras)),
                  rep(year, length(terra::cells(ras)))) %>%
    as.data.frame()
  
  names(coords) <- c("x", "y", "year") # this is *not* coord_cols
  message(names(coords))
  
  tmp <- build_design_matrix(ras,
                             coords,
                             temporal_var = FALSE,
                             scale = FALSE,
                             degs_to_rads = TRUE)
  
  X_pixel <- tmp$df %>%
    mutate(year_scaled = scaled_year)
  X_pixel <- dplyr::mutate(X_pixel, year_scaled = scaled_year)
  
  random_field_pixel <- greta.gp::project(random_field, X_pixel[,coord_cols])
  
  mut_freq_pixel <- (X_pixel[,design_cols] %*% parameters$beta + # transform here? Is R taking over?
                       random_field_pixel) %>%
    ilogit()
  
  message("here")
  
  post_pixel_sims <- greta::calculate(mut_freq_pixel,
                                      values = draws,
                                      nsim = 100,
                                      trace_batch_size = 50) # reducing: will take longer, use less mem
  message("now here")
  post_pixel_median <- apply(post_pixel_sims$mut_freq_pixel[, , 1], 2, median)
  post_pixel_sd <- apply(post_pixel_sims$mut_freq_pixel[,,1], 2, sd)
  
  out <- c(ras, ras) * 0 # assuming we have at least two layers in there ..
  names(out) <- paste0(year, c("_post_median", "_post_sd"))
  out[[1]][terra::cells(out[[1]])] <- post_pixel_median
  out[[2]][terra::cells(out[[2]])] <- post_pixel_sd
  
  if(!is.null(stable_transmission_mask)){
  #if (length(unique(suppressWarnings(values(stable_transmission_mask)))) != 1){
    # let it be known that I did some googling about this :(
    # in raster I would have plopped `raster(NA)` in the function definition
    out <- mask(out, stable_transmission_mask)
  }
  out
}



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
  
  names(coords) <- c("x", "y", "year") # this is *not* coord_cols
  
  tmp <- build_design_matrix(ras,
                             coords,
                             temporal_var = FALSE,
                             scale = FALSE,
                             degs_to_rads = TRUE)
  
  X_pixel <- tmp$df %>%
    mutate(year_scaled = scaled_year)
  
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