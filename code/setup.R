set.seed(0784)

# TODO
# filter on end.year ..
# filter on 

suppressPackageStartupMessages({
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
  library(rnaturalearthdata)})
world <- ne_countries(scale="medium", returnclass = "sf")

MIN_YEAR <- 2000
BUFFER <- 5000


nice_name_lookup <- list("k13_marcse" = "Kelch 13",
                         "crt76" = "Pfcrt-K76T",
                         "mdr86" = "Pfmdr1-N86Y",
                         "mdr184" = "Pfmdr1-Y184F",
                         "mdr1246" = "Pfmdr1-D1246Y")

data_path_lookup <- list("k13_marcse" = "data/clean/moldm_marcse_k13_nomarker.csv",
                         "crt76" = "data/clean/moldm_crt76.csv",
                         "mdr86" = "data/clean/pfmdr_single_mdr86.csv",
                         "mdr184" = "data/clean/pfmdr_single_mdr184.csv",
                         "mdr1246" = "data/clean/pfmdr_single_mdr1246.csv")


# country shps for plotting/masking
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

afr_extent <- (matrix(c(-21, -34.9, 63, 37.4), nrow = 2, ncol = 2, 
                      dimnames = list(c("x", "y"), c("min", "max"))))

# For example: 
pfpr_years <- 2000:2024
years <- as.character(pfpr_years)

# see "retrieve_pfpr.R"
pfpr <- rast("data/pfpr_rasters_afr_2025.tif") %>%
  scale()
names(pfpr) <- paste0("pfpr_", years)

covariates <- pfpr # now standardised

setup_mut_data <- function(path, min_year = NULL, buffer = 0){
  
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
    mutate(land = terra::extract(pfpr$pfpr_2000, 
                                 data.frame(x = x, y = y), 
                                 ID = FALSE,
                                 search_radius = buffer)) %>%
    drop_na() %>%
    dplyr::select(-c(land))
}

# (this was for abandoned pfmdr haplotype model:)
setup_multiple_snps <- function(path, min_year = NULL, buffer = 0){
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
    mutate(land = terra::extract(pfpr$pfpr_2000, cbind(x, y), 
                                 search_radius = buffer),
           snp_id = as.numeric(as.factor(snp))) %>%
    drop_na() %>%
    dplyr::select(-c(land))
}

# e.g.:
# mut_data <- setup_multiple_snps("../moldm/clean/pfmdr_single_locus.csv")
