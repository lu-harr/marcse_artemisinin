source("code/setup_k13.R")
source("code/build_design_matrix.R")
source("code/predict_to_raster.R")

# for prediction raster:
AGG_FACTOR = 5

scaled_years <- scale_years(range(pfpr_years))

out_dir <- "output/circmat_pfmdr86/"

# bring in all of the other outputs here too
mut_data <- read_rds(paste0(out_dir, "mut_data.rds"))
stable_transmission_mask <- rast("data/stable_transmission_mask.grd") %>%
  aggregate(AGG_FACTOR)

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

writeRaster(preds, paste0(out_dir, "preds_all.grd"), overwrite = TRUE)
