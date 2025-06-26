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


# here's how to use it:
ras <- covariates$pfpr_2000 %>%
  aggregate(fact = 2)
tmp <- survey_effort(mut_data, ras, 1.5)
writeRaster(tmp, "output/surveillance_effort.grd", overwrite = TRUE)




  
