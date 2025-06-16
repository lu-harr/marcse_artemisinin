set.seed(0784)
# for some reason, my usual seed was causing a mysterious TF error
library(greta)
#library(DiagrammeR)
library(terra)
library(readxl)
library(dplyr)
#library(malariaAtlas)
library(tidyr)
library(scales)
library(magrittr)
library(tidyverse)

afr_extent <- (matrix(c(-21, -34.9, 63, 37.4), nrow = 2, ncol = 2, 
                      dimnames = list(c("x", "y"), c("min", "max"))))

# For example: 
pfpr_years <- 2000:2022
years <- as.character(pfpr_years)

pfpr <- rast("data/pfpr_rasters_afr.tif") %>%
  scale()

names(pfpr) <- paste0("pfpr_", years)

scaled_years <- scale(pfpr_years)[,1]

scale_year <- function(year, min_year=min(pfpr_years)){
  return(scaled_years[year - min_year + 1])
}

# all_dhps_africa <- read_xls("dhfr_dhps_surveyor_data.xls",
#                             col_types=c(rep("numeric", 2), rep("guess", 2),
#                                         rep("numeric", 4), "guess",
#                                         rep("numeric", 3), "guess",
#                                         "numeric", rep("guess", 3),
#                                         "numeric", "guess", "numeric", "guess")) %>%
#   rename(x=lon, y=lat) %>%
#   filter(str_detect(`marker group`, "540E")) %>%
#   suppressWarnings()

#all_k13_africa <- read.csv("data/moldm_k13_afr_group_coord_year_tested_nomarker.csv") %>%
all_k13_africa <- read.csv("data/moldm_k13_nomarker.csv") %>% 
  rename(x = Longitude,
         y = Latitude,
         present = Present,
         tested = Tested)
  # rename(x=lon, y=lat)

zoom_k13 <- all_k13_africa %>%
  filter(x < afr_extent["x","max"] & afr_extent["x","min"] < x &
           y < afr_extent["y","max"] & afr_extent["y","min"] < y)

# {plot(pfpr$pfpr_2000, main=paste0(nrow(zoom_k13), " records in area of interest"))
#   points(as.data.frame(zoom_k13[,c("x","y")]), col="grey90",
#          cex=zoom_k13$present/zoom_k13$tested*2)
#   points(zoom_k13[zoom_k13$present/zoom_k13$tested < 0.1, c("x","y")], col="red")}

mut_data <- zoom_k13 %>%
  dplyr::select(x, y, year, tested, present) %>%
  group_by(x, y, year) %>%
  summarise(tested = sum(tested),
            present = sum(present)) %>%
  ungroup() %>%
  mutate(trunc_year = case_when(year < min(pfpr_years) ~ min(pfpr_years),
                                year > max(pfpr_years) ~ max(pfpr_years),
                                .default = year),
         scaled_year = scale_year(trunc_year, min_year=min(pfpr_years)),
         land = terra::extract(pfpr$pfpr_2000, cbind(x, y))) %>%
  drop_na() %>%
  dplyr::select(-c(land))

dim(mut_data)
# plot(mut_data$year, mut_data$present/mut_data$tested)

# {plot(pfpr$pfpr_2000, main=paste0(nrow(mut_data), " records in area of interest after grouping"))
#   points(mut_data[,c("x","y")], col="grey90",
#          cex=mut_data$present/mut_data$tested*2)
#   points(mut_data[mut_data$present/mut_data$tested < 0.1, c("x","y")], col="red")
# }

# pfpr_agg <- aggregate(pfpr, 10)
# 
# df <- cbind(xyFromCell(pfpr_agg, 1:ncell(pfpr_agg)),
#             as.matrix(pfpr_agg)) %>% # has to be a matrix
#   as.data.frame() %>%
#   pivot_longer(cols = starts_with("pfpr"), 
#                names_to = "trunc_year",
#                names_prefix = "pfpr_")

# library(rasterVis)
# ggplot(df, aes(x, y, value)) +
#   geom_raster(aes(fill = value)) +
#   geom_point(data = mut_data, 
#              colour="white", 
#              shape = 1,
#              aes(size = present/tested/2)) +
#   facet_wrap(~ trunc_year) +
#   scale_fill_gradientn(colours = hcl.colors(225)) +
#   coord_equal() +
#   labs(value = "PfPR", size="prevalence")
# # they call me a ggplot pro

covariates <- pfpr # now standardised


##########################################################
# grid_raster <- pfpr_agg[[1]] * 0
# names(grid_raster) <- NULL
# 
# # Pick projection as UTM zone 36S (Uganda, Kenya, Tanzania)
# new_crs <- "epsg:21036"
# 
# # project this template raster, the bioclim layers, and the population raster
# grid_raster <- terra::project(grid_raster, new_crs)







