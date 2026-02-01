library(malariaAtlas)
malariaAtlas::listRaster()
#  Malaria__202206_Global_Pf_Parasite_Rate 
#  Malaria__202406_Global_Pf_Parasite_Rate 
#  Malaria__202508_Global_Pf_Parasite_Rate
afr_extent <- (matrix(c(-21, -34.9, 63, 37.4), nrow = 2, ncol = 2, 
                      dimnames = list(c("x", "y"), c("min", "max"))))
pfpr_years <- 2000:2024
years <- as.character(pfpr_years)

# referred to as "covariates" everywhere else
# here's the code used to retrieve from MAP:
# previously "Malaria__202406_Global_Pf_Parasite_Rate"
pfpr <- malariaAtlas::getRaster("Malaria__202508_Global_Pf_Parasite_Rate",
                                extent=afr_extent,
                                year=years) %>% # which gives me a SpatRasterCollection
          suppressMessages() %>% # looks like some messages are still slipping the net ...
          as.list() %>%
          rast() %>%
          subset(seq(1, length(years)*2+1, 2))



writeRaster(pfpr, "data/pfpr_rasters_afr_2025.tif")


#################################################################################

# grabbing for time visualisation:
incid <- malariaAtlas::getRaster("Malaria__202508_Global_Pf_Incidence_Count",
                                 extent=afr_extent,
                                 year=2000:2024) %>% # which gives me a SpatRasterCollection
  suppressMessages() %>% # looks like some messages are still slipping the net ...
  as.list() %>%
  rast()

writeRaster(incid, "data/incid_rasters_afr.tif")

plot(incid$`Number of newly diagnosed Plasmodium falciparum cases, on a given year 2000-2024-2024`)
plot(log10(incid$`Number of newly diagnosed Plasmodium falciparum cases, on a given year 2000-2024-2024`))
# this is in the same order of magnitude as malaria report estimate:
sum(values(incid$`Number of newly diagnosed Plasmodium falciparum cases, on a given year 2000-2024-2024`), 
    na.rm = TRUE)
max(values(incid$`Number of newly diagnosed Plasmodium falciparum cases, on a given year 2000-2024-2024`),
    na.rm = TRUE)

# incid <- malariaAtlas::getRaster("Malaria__202406_Global_Pf_Incidence_Count",
#                                  extent=afr_extent,
#                                  year=years) %>% # which gives me a SpatRasterCollection
#   suppressMessages() %>% # looks like some messages are still slipping the net ...
#   as.list() %>%
#   rast()



