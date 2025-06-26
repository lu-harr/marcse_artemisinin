library(terra)
library(sf)

ras <- covariates$pfpr_2000 %>%
  aggregate(fact = 2)

survey_effort <- function(coords, # a df
                          mask, # a SpatRaster
                          sigma,
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
                          field = field,
                          fun = "sum",
                          background = 0)

  gf <- focalMat(ras, sigma, "Gauss")
  # need na.rm=TRUE or we lose edges
  test_dens <- focal(ras, gf, pad = TRUE, na.rm = TRUE)
  
  mask(test_dens, mask)
}

tmp <- survey_effort(mut_data, ras, 1.5)
writeRaster(tmp, "output/surveillance_effort.grd", overwrite = TRUE)


########################################################################
# this is all visualisation:
# 
# # mask
# afr <- world %>%
#   filter(continent == "Africa") %>%
#   vect() %>%
#   crop(ext(-21, 63, -35, 37))
# 
# # multipanel: time (hist: number tested, number of points)
# surveil <- xyFromCell(test_dens_masked, cell = cells(test_dens_masked)) %>%
#   as.data.frame() %>%
#   mutate(effort = unlist(extract(test_dens_masked, cells(test_dens_masked))))
# 
# p1 <- ggplot(data = mut_data %>%
#                group_by(year) %>%
#                summarise(n = n())) +
#   geom_bar(stat = "identity", aes(x = year, y = n)) +
#   ylab("Number of locations") +
#   xlab("Year") +
#   ggtitle("(a)")
# 
# p2 <- ggplot(data = mut_data %>%
#                group_by(year) %>%
#                summarise(tested = sum(tested))) +
#   geom_bar(stat = "identity", aes(x = year, y = tested)) +
#   ylab("Number of tests") +
#   xlab("Year") +
#   ggtitle("(b)")
# 
# # could just go back to totals? As emphasis is on effort?
# # p2 <- ggplot(data = mut_data %>%
# #                group_by(year) %>%
# #                summarise(present = sum(present),
# #                          absent = sum(tested) - present) %>%
# #                pivot_longer(!year, names_to = "Tests", values_to = "Count")) +
# #   geom_bar(stat = "identity", aes(x = year, y = Count, fill = Tests)) +
# #   xlab("Year")
# 
# p3 <- ggplot() +
#   geom_sf(data = st_as_sf(afr), fill = "white") + # not showing anything in the background here ...
#   geom_tile(data = surveil, aes(x, y, fill = effort)) +
#   geom_sf(data = st_as_sf(afr), colour = "white", fill = NA) +
#   scale_fill_viridis_c(na.value = NA, bquote(atop("Tests per","~100"~km^2)), trans="sqrt") +
#   xlab("Longitude") +
#   ylab("Latitude") +
#   ggtitle("(c)")
# 
# library(deeptime)
# 
# gg1 <- ggarrange2(p1, p2, layout = rbind(c(1), c(2)), draw = FALSE)
# ggarrange2(gg1, p3, widths = c(1,2))
# ggsave("figures/surveillance_effort.png", ggarrange2(gg1, p3, widths = c(1,2)))
# # It only took me 45 mins to work this out I guess
# 
# # Lucy: you were supposed to email HR and also the /output/ files you just tried to add
# # are too big remove from commit




  
