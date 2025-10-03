library(malariaAtlas)
malariaAtlas::listRaster()
#  Malaria__202206_Global_Pf_Parasite_Rate 
#  Malaria__202406_Global_Pf_Parasite_Rate 
afr_extent <- (matrix(c(-21, -34.9, 63, 37.4), nrow = 2, ncol = 2, 
                      dimnames = list(c("x", "y"), c("min", "max"))))
pfpr_years <- 2000:2022
years <- as.character(pfpr_years)

# referred to as "covariates" everywhere else
# here's the code used to retrieve from MAP:
pfpr <- malariaAtlas::getRaster("Malaria__202406_Global_Pf_Parasite_Rate",
                                extent=afr_extent,
                                year=years) %>% # which gives me a SpatRasterCollection
          suppressMessages() %>% # looks like some messages are still slipping the net ...
          as.list() %>%
          rast() %>%
          subset(seq(1, length(years)*2+1, 2))

writeRaster(pfpr, "data/pfpr_rasters_afr.tif")

# t

