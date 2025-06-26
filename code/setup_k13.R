set.seed(0784)

# for some reason, my usual seed was causing a mysterious TF error

library(greta)
library(greta.gp)
library(terra)
library(readxl)
library(dplyr)
library(tidyr)
library(scales)
library(magrittr)
library(tidyverse)

library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale="medium", returnclass = "sf")

# country shps for plotting/masking
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37))

afr_extent <- (matrix(c(-21, -34.9, 63, 37.4), nrow = 2, ncol = 2, 
                      dimnames = list(c("x", "y"), c("min", "max"))))

# For example: 
pfpr_years <- 2000:2022
years <- as.character(pfpr_years)

pfpr <- rast("data/pfpr_rasters_afr.tif") %>%
  scale()
# referred to as "covariates" everywhere else
covariates <- pfpr # now standardised

names(pfpr) <- paste0("pfpr_", years)

setup_mut_data <- function(path){
  # read in and format k13 data
  read.csv(path) %>% 
    rename(x = Longitude,
           y = Latitude,
           present = Present,
           tested = Tested) %>%
    filter(x < afr_extent["x","max"] & afr_extent["x","min"] < x &
             y < afr_extent["y","max"] & afr_extent["y","min"] < y) %>%
    dplyr::select(x, y, year, tested, present) %>%
    group_by(x, y, year) %>%
    summarise(tested = sum(tested),
              present = sum(present)) %>%
    ungroup() %>%
    mutate(land = terra::extract(pfpr$pfpr_2000, cbind(x, y))) %>%
    drop_na() %>%
    dplyr::select(-c(land))
}


