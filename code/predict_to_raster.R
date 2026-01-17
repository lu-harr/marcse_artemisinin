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
                           design_cols = c("year"),
                           coverage = FALSE,
                           test_pts_for_coverage = NULL){
  message("predicting to raster")
  cov_years <- names(stack)
  cov_years <- as.numeric(gsub("[^0-9]", "", cov_years))
  lab_year <- year
  if (year > max(cov_years)){
    year = max(cov_years)
  } else if (year < min(cov_years)){
    year = min(cov_years)
  }
  
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
  
  # temporal var is turned off ... all one year
  tmp <- build_design_matrix(ras,
                             coords,
                             temporal_var = FALSE,
                             scale = FALSE,
                             degs_to_rads = TRUE)
  
  # finish off design matrix
  X_pixel <- tmp$df %>%
    dplyr::mutate(year_scaled = scaled_year)
  #  message(paste0("Scaled year: ", scaled_year))
  # X_pixel <- dplyr::mutate(X_pixel, year_scaled = scaled_year)
  # message(coord_cols) # this is xyyear
  # message(design_cols) # this is interceptyear_scaledpfpr
  # message(paste(names(X_pixel), collapse = ", "))
  # message(dim(X_pixel[,coord_cols]))
  # message(dim(X_pixel[,design_cols]))
  # 
  # project random field to coordinates we would like predictions for
  random_field_pixel <- greta.gp::project(random_field, X_pixel[,coord_cols])
  
  # message(dim(t(parameters$beta)))
  # message(dim(random_field_pixel))
  mut_freq_pixel <- (X_pixel[,design_cols] %*% parameters$beta + 
                       # transform here? Is R taking over?
                       random_field_pixel) %>%
    ilogit()
  
  # to do: have a look at variance/lengthscale trade-off
  
  post_pixel_sims <- greta::calculate(mut_freq_pixel,
                                      values = draws,
                                      nsim = nsim,
                                      trace_batch_size = 1) # reducing: will take longer, use less mem
  
  if (coverage & !is.null(test_pts_for_coverage)){
    coverages <- calculate_coverages(post_pixel_sims,
                                      test_pts_for_coverage,
                                      lab_year,
                                      ras,
                                     incs = 100)
  } else {
    coverages <- NULL
  }
  
  # message(coverages)
  
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
  message("predicted to raster")
  list(out = out, coverages = coverages)
}

# out_dir <- "output/k13_marcse/gneiting_sparse/"
# source("code/setup.R")
# source("code/build_design_matrix.R")
# scaled_years <- scale_years(range(pfpr_years))
# 
# # bring in all of the other outputs here too
# AGG_FACTOR <- 10
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

predict_to_points <- function(pts, 
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
                               coverage = FALSE,
                               data_path = ""){
  message("predicting to points")
  message("warning coverages call is deprecated")
  
  # coords should be npts * c("x", "y", "year")
  # X_pixel should be 
  
  
  cov_years <- names(stack)
  cov_years <- as.numeric(gsub("[^0-9]", "", cov_years))
  lab_year <- year
  if (year > max(cov_years)){
    year = max(cov_years)
  } else if (year < min(cov_years)){
    year = min(cov_years)
  }
  
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
  
  # temporal var is turned off ... all one year
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
  
  # message(ncell(ras))
  
  post_pixel_sims <- greta::calculate(mut_freq_pixel,
                                      values = draws,
                                      nsim = nsim,
                                      trace_batch_size = 1) # reducing: will take longer, use less mem
  
  if (coverage & data_path != ""){
    coverages <- calculate_coverages(post_pixel_sims,
                                     data_path,
                                     lab_year,
                                     ras,
                                     incs = 100)
  } else {
    coverages <- NULL
  }
  
  # message(coverages)
  
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
  message("predicted to raster")
  list(out = out, coverages = coverages)
}






#' Calculate coverage probabilities
#'
#' @param sims predictions from greta::calculate()
#' @param path character. Path to mut_data object
#' @param yr numeric. Year predictions are relevant to.
#' @param ras SpatRaster. Land mask for predictions.
#' @param incs numeric. Number of increments in credible interval widths to make in [0,100]
#'
#' @returns list of length incs + 3: incs + entries for coverage probabilities; 
#' `n` for number of entries in `mut_data` for `yr`; `n_landed` for number of entries (of `n`)
#' that are at non-NA values of `ras`.
#' @export
#'
#' @examples
calculate_coverages <- function(sims, dat, yr, ras, incs = 100){
  out <- list()
  
  dat <- dat %>%
    filter(year == yr)
  
  if (nrow(dat) == 0){
    return(NULL)
  }
    
  # now get indices for matrix
  dat <- dat %>%
    mutate(cell = terra::extract(ras, dat[,c("x", "y")], cells = TRUE) %>%
                  dplyr::select(cell) %>%
                  unlist()) %>%
    left_join(data.frame(cell = cells(ras), 
                         idx = 1:length(cells(ras))),
                   by = join_by(cell == cell))
  
  out$n_pts <- nrow(dat)
  out$n_landed <- sum(!is.na(dat$idx))
  out$n_non_zero <- sum(!is.na(dat$idx) & dat$present > 0)
  
  # could add a re-landing step
  dat <- filter(dat, !is.na(idx))
  
  widths = seq(0, 100, length.out = incs + 1)
  probs = unique(c(0.5 - widths / 200, 0.5 + widths / 200))
  
  idx <- unique(dat$idx)
  
  # message(length(idx))
  # message(paste("dim", dim(sims$mut_freq_pixel)))

  if (length(idx) == 0){
    message("no coverages to check")
    return(NULL)
  }


  if (length(idx) == 1){
    message("Here's a fix")
    bounds <- quantile(sims$mut_freq_pixel[,idx,1], probs = probs) %>%
      t() %>%
      as.data.frame() %>%
      mutate(idx = idx)
  } else {
    # add cell IDs into middle index here and save yourself some time
    bounds <- apply(sims$mut_freq_pixel[,idx,1], 2, quantile,
                      probs = probs) %>% # a nprobs * ncell matrix
        t() %>%
        as.data.frame() %>%
        mutate(idx = idx)
  }
  
  dat <- left_join(dat, bounds, by = join_by(idx)) %>%
    mutate(prevalence = present / tested)
  
  dat_non_zero <- dat %>%
    filter(prevalence > 0)
  
  # message(paste(names(dat)))
  
  ind <<- 9 # number of columns to skip
  coverages <- lapply(widths, function(width){
    # pick out the upper and lower bound columns for this width ..
    # this probably a nice way to do this based on column ordering ..
    if (width == 0){return(0)}
    #lower <- which(names(dat) == paste0(as.character(50 - width / 2), "%"))
    #upper <- which(names(dat) == paste0(as.character(50 + width / 2), "%"))
    lower <- ind
    upper <- ind + incs
    out <- c(sum(dat$prevalence >= dat[,lower] & dat$prevalence <= dat[,upper]),
             sum(dat_non_zero$prevalence >= dat[,lower] & dat_non_zero$prevalence <= dat[,upper]))
    # message(paste(width, lower, upper, out, ind))
    ind <<- ind + 1
    out
  })
  names(coverages) <- widths
  
  c(out, coverages)
}

# as in predict_to_ras:
# tmp = calculate_coverages(post_pixel_sims,
#                     "output/k13_marcse/gneiting_sparse/",
#                     2022,
#                     ras, incs = 100)


#' Concatenate multiple (annual) coverage sheets together
#'
#' @param path character
#' @param fold numeric - for cross validation
#'
#' @returns
#' @export
#'
#' @examples
concat_coverages <- function(path, fold = NULL){
  
  if(is.null(fold)){
    path <- paste0(path, "coverages/")
    fold <- ""
  }

  files <- path %>% list.files()
  
  if(length(grep(pattern = "\\d{4}_coverages", files)) == 0){
    message("Nothing to concat!")
    return(0)
  }
  files <- paste0(path, files[grepl(pattern = "\\d{4}_coverages", files)])

  if(fold != ""){
    files <- files[str_extract(files, ".(?=\\.csv$)") == as.character(fold)]
  }

  # open all of the files and rbind them together
  coverages <- sapply(files, read.csv) %>%
    rbind() %>%
    t() %>%
    as.data.frame() %>%
    # retrieve year from rownames
    mutate(year = str_extract(rownames(.), "\\d{4}")) %>%
    unnest(cols = everything())
  
  write.csv(coverages, paste0(path, "all_coverages", fold, ".csv"), row.names = FALSE)
  message("coverages concatted")
  
  return(str_extract(files, "\\d{4}"))
}

# use cases:
# concat_coverages("output/k13_marcse/bb_gne/")
# 
# concat_coverages("output/k13_marcse/bb_gne/cv_preds/", 9)

# write annual preds objects as medians/sds preds objects
concat_preds <- function(path, 
                         medians = TRUE, 
                         sds = FALSE, 
                         sdscaled = FALSE, 
                         ciwidth = FALSE,
                         upper = FALSE, 
                         lower = FALSE,
                         fold = ""){
  if (fold != ""){fold <- paste0("_", fold)}
  
  to_read <- grep(paste0("\\d{4}_preds", fold, ".grd$"), 
                  list.files(path), value = TRUE)
  
  if(length(to_read) == 0){
    message("Nothing to concat!")
    return(0)
  }
  
  razzes <- rast(paste0(path, "/", to_read))
  
  to_write <- c()

  if (medians){
    to_write <- names(razzes)[grep("50", names(razzes))]
    
    terra::writeRaster(subset(razzes, to_write), 
                        paste0(path, "preds_medians", fold, ".tif"), 
                        overwrite = TRUE, filetype = "GTiff")
  }
  if (sds){
    to_write <- names(razzes)[grep("sd$", names(razzes))]
    terra::writeRaster(subset(razzes, to_write), 
                        paste0(path, "preds_sds", fold, ".tif"), 
                        overwrite = TRUE, filetype = "GTiff")
  }
  if (sdscaled){
    to_write <- names(razzes)[grep("sdscaled", names(razzes))]
    terra::writeRaster(subset(razzes, to_write), 
                        paste0(path, "preds_sdscaled", fold, ".tif"), 
                        overwrite = TRUE, filetype = "GTiff")
  }
  if (upper){
    to_write <- names(razzes)[grep("97\\.5", names(razzes))]
    terra::writeRaster(subset(razzes, to_write), 
                       paste0(path, "preds_upper", fold, ".tif"), 
                        overwrite = TRUE, filetype = "GTiff")
  }
  if (lower){
    to_write <- names(razzes)[grep("2\\.5", names(razzes))]
    terra::writeRaster(subset(razzes, to_write), 
                       paste0(path, "preds_lower", fold, ".tif"), 
                        overwrite = TRUE, filetype = "GTiff")
  }
  if (ciwidth){
    ci_width = subset(razzes, grepl("97\\.5", names(razzes))) - subset(razzes, grepl("2\\.5", names(razzes)))
    names(ci_width) = paste0(str_extract(names(ci_width), "\\d{4}"), "_CI")
    terra::writeRaster(ci_width, 
                       paste0(path, "preds_ciwidths", fold, ".tif"), 
                        overwrite = TRUE, filetype = "GTiff")
  }

  # message("\n")
  # message(paste(to_write))
  
  # razzes <- subset(razzes, to_write)
  # f <- file.path(path, "preds_all.tif")
  # terra::writeRaster(razzes, f, overwrite = TRUE, filetype = "GTiff")
  message("preds concatted")
}

# use cases
# concat_preds("output/k13_marcse/bb_gne/",
#              medians = TRUE)
# 
# concat_preds("output/k13_marcse/bb_gne/cv_preds/", 
#              medians = TRUE, fold = 1)




