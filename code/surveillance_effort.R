survey_effort <- function(coords, # a df
                          sigma,
                          mask, # a SpatRaster
                          actually_really_apply_mask = TRUE,
                          xy = c("x", "y"), # where the coords are
                          field = "tested", # where the counts are
                          crs = "epsg:4326"){
  # given df of `coords` associated with column `field` (by default, number of 
  # tests), find Gaussian kernel density estimate with standard deviation `sigma`
  # return SpatRaster masked by `mask` (i.e. `stable_transmission_mask`)
  
  # cast to SpatVector
  coords <- coords %>%
    vect(geom = xy, crs = crs) %>%
    terra::project(mask)
  
  # sum "tested" in each pixel
  ras <- terra::rasterize(coords,
                          mask,
                          field = field, # sum of number of tests in pixel
                          fun = "sum",
                          background = 0)

  gf <- focalMat(ras, sigma, "Gauss")
  # need na.rm=TRUE or we lose edges
  test_dens <- focal(ras, gf, pad = TRUE, na.rm = TRUE)
  
  if (!is.null(mask)){test_dens <- mask(test_dens, mask)}
  
  test_dens
}



wrap_survey_effort <- function(mut_data_path,
                               out_path,
                               mask_path = "data/stable_transmission_mask.grd",
                               apply_mask = FALSE,
                               agg_factor = 1,
                               sigma = 1.5,
                               years = NULL,
                               bin_years = FALSE,
                               plot_out = NULL){
  print(out_path)
  # not 100% sure why I'm wrapping a wrapper ...
  # A wrapper for `survey_effort()` that reads in `coords` and `mask`, does optional
  # aggregation of `mask`, optionally separates `coords` by year and makes multiple
  # (annual) calls to `survey_effort()`.
  if (!is.null(mask_path)){
    transmission_mask <- rast(mask_path)
    if (agg_factor > 1){
      transmission_mask <- aggregate(transmission_mask, agg_factor)
    }
  }
  
  mut_data <- setup_mut_data(mut_data_path)
  
  
  if (!is.null(years)){
    test_dens <- lapply(1:(length(years) - as.numeric(bin_years)), function(i){
      
      # tried to do this with an ifelse which was a demoralising waste of time
      if (bin_years) {
        tmp <- filter(mut_data, year >= years[i] & year < years[i + 1])
      } else {
        tmp <- mut_data %>% filter(year == years[i])
      }
      
      survey_effort(coords = tmp, 
                    mask = transmission_mask,
                    actually_really_apply_mask = apply_mask,
                    sigma = sigma)
      
    })
    test_dens <- rast(test_dens)
    names(test_dens) <- paste0("surv_", years[1:(length(years) - as.numeric(bin_years))])
  } else {
    test_dens <- survey_effort(mut_data, 
                               mask = transmission_mask,
                               actually_really_apply_mask = apply_mask,
                               sigma = sigma)
  }
  
  if (!is.null(plot_out)){plot(test_dens)}
  
  writeRaster(test_dens, out_path, overwrite = TRUE)
}

tmp = wrap_survey_effort("data/clean/moldm_k13_nomarker.csv",
                   "output/k13/surveillance_effort_k13.grd",
                   years = seq(2012, 2024, 3), # watch out for 2025
                   bin_years = TRUE,
                   sigma = 1.5, apply_mask = FALSE)

###############################################################################
library(terra)
library(sf)
source("code/setup.R")

# function calls:
wrap_survey_effort("data/clean/moldm_k13_nomarker.csv",
                   "output/k13/surveillance_effort_k13.grd",
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/moldm_k13_nomarker.csv",
                   "output/k13/surveillance_effort_k13_tempo.grd",
                   years = seq(2012, 2024, 3),
                   bin_years = TRUE,
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/moldm_crt76.csv",
                   "output/crt76/surveillance_effort_crt.grd",
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/moldm_crt76.csv",
                   "output/crt76/surveillance_effort_crt_tempo.grd",
                   years = seq(2012, 2024, 3),
                   bin_years = TRUE,
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/pfmdr_single_86.csv",
                   "output/mdr86/surveillance_effort_pfmdr86.grd",
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/pfmdr_single_184.csv",
                   "output/mdr184/surveillance_effort_pfmdr184.grd",
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/pfmdr_single_1246.csv",
                   "output/mdr1246/surveillance_effort_pfmdr1246.grd",
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/pfmdr_single_86.csv",
                   "output/mdr86/surveillance_effort_pfmdr86_tempo.grd",
                   years = seq(2012, 2024, 3),
                   bin_years = TRUE,
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/pfmdr_single_184.csv",
                   "output/mdr184/surveillance_effort_pfmdr184_tempo.grd",
                   years = seq(2012, 2024, 3),
                   bin_years = TRUE,
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/pfmdr_single_1246.csv",
                   "output/mdr1246/surveillance_effort_pfmdr1246_tempo.grd",
                   years = seq(2012, 2024, 3),
                   bin_years = TRUE,
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/pfmdr_single_locus.csv",
                   "output/mdr_hier/surveillance_effort_pfmdragg_tempo.grd",
                   years = seq(2012, 2024, 3),
                   bin_years = TRUE,
                   sigma = 1.5, apply_mask = FALSE)

# now with added unpublished MARCSE data from dashboard paper:
wrap_survey_effort("data/clean/moldm_marcse_k13_nomarker.csv",
                   "output/k13_marcse/surveillance_effort_k13.grd",
                   sigma = 1.5, apply_mask = FALSE)

wrap_survey_effort("data/clean/moldm_marcse_k13_nomarker.csv",
                   "output/k13_marcse/surveillance_effort_k13_tempo.grd",
                   years = seq(2012, 2024, 3),
                   bin_years = TRUE,
                   sigma = 1.5, apply_mask = FALSE)
