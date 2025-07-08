source("code/setup.R")
source("code/build_design_matrix.R")
source("code/predict_to_raster.R")

# for prediction raster:
AGG_FACTOR = 5

# set this to the location where all the inference outputs are at:
out_dir <- "output/circmat_k13/"

scaled_years <- scale_years(range(pfpr_years))

# bring in all of the other outputs here too
mut_data <- read_rds(paste0(out_dir, "mut_data.rds"))
stable_transmission_mask <- rast("data/stable_transmission_mask.grd") %>%
  aggregate(AGG_FACTOR)
random_field <- read_rds(paste0(out_dir, "random_field.rds"))
parameters <- read_rds(paste0(out_dir, "parameters.rds"))
draws <- read_rds(paste0(out_dir, "draws.rds"))

pfpr_years = 2000:2022

set.seed(0748)
preds <- lapply(pfpr_years, function(year){
  message(year)
  predict_to_ras(covariates,
                 year,
                 draws,
                 parameters,
                 random_field,
                 agg_factor = AGG_FACTOR,
                 scaled_year = scaled_years[[as.character(year)]],
                 coord_cols = c("x_rd", "y_rd", "year_scaled"),
                 design_cols = c("intercept","year_scaled","pfpr"),
                 stable_transmission_mask = stable_transmission_mask)
})

preds <- rast(preds)

writeRaster(preds, paste0(out_dir, "preds_all.grd"), overwrite = TRUE)
