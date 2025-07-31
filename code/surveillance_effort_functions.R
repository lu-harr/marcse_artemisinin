#' Calculate Gaussian KDE for data.frame of surveillance data: numbers of tests 
#' associated with lons and lats
#'
#' @param coords a data.frame
#' @param sigma numeric: bandwidth/standard deviation for KDE
#' @param mask SpatRaster: e.g. stable_transmission_mask. Applied post-hoc
#' @param actually_really_apply_mask boolean: actually really apply mask? just me being silly
#' @param xy vector of chars: names of columns in `coords` corresponding to lons and lats
#' @param field character: name of column in `coords` where survey effort is, e.g. `"tested"`
#' @param crs character
#'
#' @returns SpatRaster
#' @export
#'
#' @examples
survey_effort <- function(coords, # a df
                          sigma,
                          mask, # a SpatRaster
                          actually_really_apply_mask = TRUE,
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
                          field = field, # sum of number of tests in pixel
                          fun = "sum",
                          background = 0)

  gf <- focalMat(ras, sigma, "Gauss")
  
  # need na.rm=TRUE or we lose edges
  test_dens <- focal(ras, gf, pad = TRUE, na.rm = TRUE)
  
  if (!is.null(mask) & actually_really_apply_mask){test_dens <- mask(test_dens, mask)}
  
  test_dens
}


#' Wrap calls to `survey_effort()`, e.g. for year binning
#'
#' @param mut_data_path path to mut_data
#' @param out_path path to where we should put output raster/s
#' @param mask_path path to mask
#' @param apply_mask boolean: apply mask?
#' @param agg_factor numeric: option to aggregate mask for lower resolution KDE
#' @param sigma standard deviation of Gaussian KDE/s
#' @param years vector of numerics; if `bin_years` is TRUE, then should include lower 
#' and upper bounds (low inclusive, upper non-inclusive)
#' @param bin_years boolean: if TRUE, then bin data by year, using years in `years` 
#' as breaks; if FALSE, then only calculate KDE/s for specific years in `years`
#' @param plot_out boolean: plot output ?
#' @param msg_me boolean: message year to std_out before each call to `survey_effort()`
#'
#' @returns
#' @export
#'
#' @examples
wrap_survey_effort <- function(mut_data_path,
                               out_path = NULL,
                               mask_path = "data/stable_transmission_mask.grd",
                               apply_mask = FALSE,
                               agg_factor = 1,
                               sigma = 1.5,
                               years = NULL,
                               bin_years = FALSE,
                               plot_out = FALSE,
                               msg_me = FALSE){
  
  if (!is.null(mask_path)){
    transmission_mask <- rast(mask_path)
    if (agg_factor > 1){
      transmission_mask <- aggregate(transmission_mask, agg_factor)
    }
  }
  
  mut_data <- setup_mut_data(mut_data_path)
  
  
  if (!is.null(years)){
    test_dens <- lapply(1:(length(years) - as.numeric(bin_years)), function(i){
      # the as.numeric(bin_years) skips the last item in the vector when we're 
      # binning *between* years
      # tried to do this with an ifelse which was a demoralising waste of time
      if (bin_years) {
        tmp <- filter(mut_data, year >= years[i] & year < years[i + 1])
      } else {
        tmp <- filter(mut_data, year == years[i])
      }
      
      if (msg_me) {message(years[i])}
      
      survey_effort(coords = tmp, 
                    mask = transmission_mask,
                    actually_really_apply_mask = apply_mask,
                    sigma = sigma)
      
    })
    test_dens <- rast(test_dens)
    
    # let's make this nice
    if (!bin_years){
      names(test_dens) <- paste0("surv_", years)
    } else {
      names(test_dens) <- paste0("surv_", years[1:(length(years) - 1)], "_", years[2:length(years)])
    }
    
  } else {
    test_dens <- survey_effort(mut_data, 
                               mask = transmission_mask,
                               actually_really_apply_mask = apply_mask,
                               sigma = sigma)
  }
  
  if (plot_out){plot(test_dens)}
  
  if (!is.null(out_path)){writeRaster(test_dens, out_path, overwrite = TRUE)}
  
  test_dens
}

# e.g.:
# tmp = wrap_survey_effort("data/clean/moldm_k13_nomarker.csv",
#                          sigma = 1.5, apply_mask = FALSE, agg_factor = 10)
# 
# tmp = wrap_survey_effort("data/clean/moldm_k13_nomarker.csv",
#                          years = seq(2012, 2024, 3),
#                          agg_factor = 10,
#                          bin_years = TRUE,
#                          sigma = 1.5, apply_mask = FALSE)
# 
# tmp = wrap_survey_effort("data/clean/moldm_k13_nomarker.csv",
#                          years = seq(2012, 2024, 3),
#                          agg_factor = 10,
#                          bin_years = FALSE,
#                          sigma = 1.5, apply_mask = FALSE)

###############################################################################
# library(terra)
# library(sf)
# source("code/setup.R")
# 
# # function calls:
# wrap_survey_effort("data/clean/moldm_k13_nomarker.csv",
#                    "output/k13/surveillance_effort_k13.grd",
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/moldm_k13_nomarker.csv",
#                    "output/k13/surveillance_effort_k13_tempo.grd",
#                    years = seq(2012, 2024, 3),
#                    bin_years = TRUE,
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/moldm_crt76.csv",
#                    "output/crt76/surveillance_effort_crt.grd",
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/moldm_crt76.csv",
#                    "output/crt76/surveillance_effort_crt_tempo.grd",
#                    years = seq(2012, 2024, 3),
#                    bin_years = TRUE,
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/pfmdr_single_86.csv",
#                    "output/mdr86/surveillance_effort_pfmdr86.grd",
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/pfmdr_single_184.csv",
#                    "output/mdr184/surveillance_effort_pfmdr184.grd",
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/pfmdr_single_1246.csv",
#                    "output/mdr1246/surveillance_effort_pfmdr1246.grd",
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/pfmdr_single_86.csv",
#                    "output/mdr86/surveillance_effort_pfmdr86_tempo.grd",
#                    years = seq(2012, 2024, 3),
#                    bin_years = TRUE,
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/pfmdr_single_184.csv",
#                    "output/mdr184/surveillance_effort_pfmdr184_tempo.grd",
#                    years = seq(2012, 2024, 3),
#                    bin_years = TRUE,
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/pfmdr_single_1246.csv",
#                    "output/mdr1246/surveillance_effort_pfmdr1246_tempo.grd",
#                    years = seq(2012, 2024, 3),
#                    bin_years = TRUE,
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/pfmdr_single_locus.csv",
#                    "output/mdr_hier/surveillance_effort_pfmdragg_tempo.grd",
#                    years = seq(2012, 2024, 3),
#                    bin_years = TRUE,
#                    sigma = 1.5, apply_mask = FALSE)
# 
# # now with added unpublished MARCSE data from dashboard paper:
# wrap_survey_effort("data/clean/moldm_marcse_k13_nomarker.csv",
#                    "output/k13_marcse/surveillance_effort_k13.grd",
#                    sigma = 1.5, apply_mask = FALSE)
# 
# wrap_survey_effort("data/clean/moldm_marcse_k13_nomarker.csv",
#                    "output/k13_marcse/surveillance_effort_k13_tempo.grd",
#                    years = seq(2012, 2024, 3),
#                    bin_years = TRUE,
#                    sigma = 1.5, apply_mask = FALSE)
