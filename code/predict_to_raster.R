# why no roxygen?

predict_to_ras <- function(stack, 
                           year, 
                           draws, 
                           parameters,
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
  
  names(coords) <- c("x", "y", "year")
  
  tmp <- build_design_matrix(ras,
                             coords,
                             temporal_var = FALSE,
                             scale = FALSE,
                             degs_to_rads = TRUE)
  
  X_pixel <- tmp$df %>%
    mutate(year_scaled = scaled_year)
  
  # Rand field is coming from above
  random_field_pixel <- greta.gp::project(random_field, X_pixel[,coord_cols])
  
  mut_freq_pixel <- (X_pixel[,design_cols] %*% parameters$beta + random_field_pixel) %>%
    ilogit()
  
  post_pixel_sims <- greta::calculate(mut_freq_pixel,
                                      values = draws,
                                      nsim = 100,
                                      trace_batch_size = 50) # reducing: will take longer, use less mem

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

# predict_to_ras(covariates,
#                2010,
#                draws,
#                parameters,
#                agg_factor = 15,
#                stable_transmission_mask = stable_transmission_mask)
