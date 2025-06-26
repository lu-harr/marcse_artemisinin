#source("code/setup.R")
# need to retrieve scaled_years FROM SOMEWHERE THAT IS NOT SETUP SCRIPT
source("code/build_design_matrix.R")
source("code/predict_to_raster.R")

# bring in all of the other outputs here too
mut_data <- read_rds("output/circmat_model/mut_data.rds")
stable_transmission_mask <- rast("data/stable_transmission_mask.grd")

preds <- lapply(pfpr_years, function(year){
  message(year)
  predict_to_ras(covariates,
                 year,
                 draws,
                 parameters,
                 agg_factor = AGG_FACTOR,
                 scaled_year = scaled_years[[as.character(year)]],
                 coord_cols = c("x_rd", "y_rd", "year_scaled"),
                 design_cols = c("intercept","year_scaled","pfpr"),
                 stable_transmission_mask = stable_transmission_mask)
})

preds <- rast(preds)
# made edit to pred_to_ras
# names(preds) <- sort(apply(expand.grid(pfpr_years, c("post_mean", "post_sd")), 1, 
#                          paste, collapse = "_"))

writeRaster(preds, "output/circmat_model/preds_all.grd", overwrite = TRUE)
