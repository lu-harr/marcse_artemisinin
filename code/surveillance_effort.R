library(terra)
library(sf)

survey_effort <- function(coords, # a df
                          mask, # a SpatRaster
                          sigma,
                          xy = c("x", "y"), # where the coords are
                          field = "tested", # where the counts are
                          crs = "epsg:4326"){
  # cast to SpatVector
  coords <- coords %>%
    vect(geom = xy, crs = crs) %>%
    terra::project(mask)
  
  # sum "tested" in each pixel
  ras <- terra::rasterize(coords,
                          mask,
                          field = field,
                          fun = "sum",
                          background = 0)

  gf <- focalMat(ras, sigma, "Gauss")
  # need na.rm=TRUE or we lose edges
  test_dens <- focal(ras, gf, pad = TRUE, na.rm = TRUE)
  
  mask(test_dens, mask)
}

wrap_survey_effort <- function(mut_data_path,
                               out_path,
                               agg_factor = 1,
                               sigma = 1.5,
                               years = NULL,
                               plot_out = NULL){
  # not 100% sure why I'm wrapping a wrapper
  stable_transmission_mask <- rast("data/stable_transmission_mask.grd")
  if (agg_factor > 1){
    stable_transmission_mask <- aggregate(stable_transmission_mask, agg_factor)
  }
  
  mut_data <- setup_mut_data(mut_data_path)
  
  if (!is.null(years)){
    test_dens <- lapply(years, function(z){
      tmp <- mut_data %>% filter(year == z)
      survey_effort(tmp, stable_transmission_mask, 1.5)
    })
    test_dens <- rast(test_dens)
    names(test_dens) <- paste0("surv_", years)
  } else {
    test_dens <- survey_effort(mut_data, stable_transmission_mask, sigma)
  }
  
  if (!is.null(plot_out)){plot(test_dens)}
  writeRaster(test_dens, out_path, overwrite = TRUE)
}


wrap_survey_effort("data/moldm_k13_nomarker.csv",
                   "output/circmat_k13/surveillance_effort_k13.grd",
                   sigma = 1.5, 
                   plot_out = TRUE)

wrap_survey_effort("data/moldm_k13_nomarker.csv",
                   "output/circmat_k13/surveillance_effort_k13_testyr.grd",
                   years = 2020:2021,
                   agg_factor = 2,
                   sigma = 1.5, 
                   plot_out = TRUE)

#writeRaster(test_dens, "output/circmat_k13/surveillance_effort_k13_2010-23.grd", overwrite = TRUE)

# For pfcrt:
mut_data <- setup_mut_data("data/moldm_crt76.csv")
test_dens <- survey_effort(mut_data, stable_transmission_mask, 1.5)
plot(test_dens, main = "Pfcrt76 surveillance effort")
writeRaster(test_dens, "output/circmat_crt/surveillance_effort_crt.grd", overwrite = TRUE)

# by year:
years = 2010:2023
test_dens <- lapply(years, function(z){
  tmp <- mut_data %>% filter(year == z)
  survey_effort(tmp, stable_transmission_mask, 1.5)
})

test_dens <- rast(test_dens)
names(test_dens) <- paste0("surv_", years)
plot(test_dens)
writeRaster(test_dens, "output/circmat_crt/surveillance_effort_crt_2010-23.grd", overwrite = TRUE)




  
