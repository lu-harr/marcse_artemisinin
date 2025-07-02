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

stable_transmission_mask <- rast("data/stable_transmission_mask.grd")

# here's how to use it:
# ras <- covariates$pfpr_2000 %>%
#   aggregate(fact = 2)

# For K13:
mut_data <- setup_mut_data("data/moldm_k13_nomarker.csv")
test_dens <- survey_effort(mut_data, stable_transmission_mask, 1.5)
plot(test_dens, main = "Kelch 13 surveillance effort")
writeRaster(test_dens, "output/circmat_k13/surveillance_effort_k13.grd", overwrite = TRUE)

# by year:
years = 2010:2023
test_dens <- lapply(years, function(z){
  tmp <- mut_data %>% filter(year == z)
  survey_effort(tmp, stable_transmission_mask, 1.5)
})

test_dens <- rast(test_dens)
names(test_dens) <- paste0("surv_", years)
plot(test_dens)
writeRaster(test_dens, "output/circmat_k13/surveillance_effort_k13_2010-23.grd", overwrite = TRUE)

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




  
