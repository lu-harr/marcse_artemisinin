set.seed(0784)

# TODO
# filter on end.year ..
# filter on 

library(greta)
library(greta.gp)
library(terra)
library(readxl)
library(dplyr)
library(tidyr)
library(scales)
library(magrittr)
library(tidyverse)
library(sf) # not sure why I'm adding this in October 2025
library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale="medium", returnclass = "sf")

MIN_YEAR <- 2000

# country shps for plotting/masking
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

afr_extent <- (matrix(c(-21, -34.9, 63, 37.4), nrow = 2, ncol = 2, 
                      dimnames = list(c("x", "y"), c("min", "max"))))

# For example: 
pfpr_years <- 2000:2022
years <- as.character(pfpr_years)

# see "retrieve_pfpr.R"
pfpr <- rast("data/pfpr_rasters_afr.tif") %>%
  scale()
names(pfpr) <- paste0("pfpr_", years)

covariates <- pfpr # now standardised

setup_mut_data <- function(path, min_year = NULL){
  
  # read in and format data for a single response (k13, crt76, mdr1-86, ...)
  read.csv(path) %>% 
    rename(x = Longitude,
           y = Latitude,
           present = Present,
           tested = Tested) %>%
    filter(if (is.null(min_year)) TRUE else year >= min_year) %>%
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

setup_multiple_snps <- function(path, min_year = NULL){
  # read in and format data for multiple responses
  # (extra columns: snp, snp_id; snp included in grouping)
  read.csv(path) %>% 
    rename(x = Longitude,
           y = Latitude,
           present = Present,
           tested = Tested,
           snp = realised) %>%
    # remove floaters outside of rectangular extent
    filter(x < afr_extent["x","max"] & afr_extent["x","min"] < x &
             y < afr_extent["y","max"] & afr_extent["y","min"] < y) %>%
    dplyr::select(x, y, year, tested, present, snp) %>%
    group_by(x, y, year, snp) %>%
    summarise(tested = sum(tested),
              present = sum(present)) %>%
    ungroup() %>%
    # check all points are landed
    mutate(land = terra::extract(pfpr$pfpr_2000, cbind(x, y)),
           snp_id = as.numeric(as.factor(snp))) %>%
    drop_na() %>%
    dplyr::select(-c(land))
}

# e.g.:
# mut_data <- setup_multiple_snps("../moldm/clean/pfmdr_single_locus.csv")
