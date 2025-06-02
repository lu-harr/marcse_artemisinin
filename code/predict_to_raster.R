predict_to_ras <- function(stack, 
                           year, 
                           draws, 
                           parameters,
                           agg_factor = 1,
                           nsim = 100,
                           stable_transmission_mask = rast()){
  ras <- stack[[grep(as.character(year), names(stack))]] %>%
    aggregate(fact = agg_factor)
  
  
  coords_pixel <- cbind(terra::xyFromCell(ras, cell = terra::cells(ras)),
                        rep(year, length(terra::cells(ras)))) %>%
    as.data.frame() %>%
    setNames(c("x", "y", "year")) %>%
    mutate(scaled_year = scale_year(year))
  message(nrow(coords_pixel))
  
  X_pixel <- build_design_matrix(stack,
                                 coords_pixel,
                                 temporal_range = rep(year, 2),
                                 scale = FALSE)
  message(nrow(X_pixel))
  
  coords_pixel <- dplyr::select(coords_pixel, -c("year"))
  
  # Rand field is coming from above
  random_field_pixel <- greta.gp::project(random_field, coords_pixel)
  
  mut_freq_pixel <- (X_pixel %*% parameters$beta + random_field_pixel) %>%
    ilogit()
  
  message(dim(mut_freq_pixel))
  
  post_pixel_sims <- greta::calculate(mut_freq_pixel,
                                      values = draws,
                                      nsim = 100,
                                      trace_batch_size = 10) # reducing: will take longer, use less mem
  message(dim(mut_freq_pixel))
  
  post_pixel_mean <- colMeans(post_pixel_sims$mut_freq_pixel[, , 1])
  post_pixel_sd <- apply(post_pixel_sims$mut_freq_pixel[,,1], 2, sd)
  
  out <- c(ras, ras) * 0 # assuming we have at least two layers in there ..
  out <- setNames(out, c("post_mean", "post_sd"))
  
  out$post_mean[terra::cells(out$post_mean)] <- post_pixel_mean
  out$post_sd[terra::cells(out$post_sd)] <- post_pixel_sd
  
  if (length(unique(suppressWarnings(values(stable_transmission_mask)))) != 1){
    # let it be known that I did some googling about this :(
    # in raster I would have plopped `raster(NA)` in the function definition
    out <- mask(out, stable_transmission_mask)
  }
  
  out
}

# predict_to_ras(covariates,
#                2010,
#                draws,
#                parameters,
#                agg_factor = 15,
#                stable_transmission_mask = stable_transmission_mask)
