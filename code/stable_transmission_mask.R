# cook a stable transmission mask ... 
# subject to change but would like to be smoothish
# (get feedback from Jen/Philippe/Karen before publication obvs)
source("code/setup.R")
library(sf)

ras <- covariates$pfpr_2021 %>%
  aggregate(AGG_FACTOR)

wt <- focalMat(ras, 0.5, type = "rectangle")
ras_blur <- focal(ras - minmax(ras)[1], wt, na.rm=TRUE)

ras_mask <- ras_blur
ras_mask[ras_mask < 0.1] <- NA
ras_mask <- mask(ras_mask, ras)
ras_mask <- mask(ras_mask, 
                 world %>% 
                   filter(continent == "Africa") %>% 
                   st_buffer(dist = 100000))

# plot(ras_blur)
# plot(ras, col="grey80")
# plot(ras_mask, add=TRUE)
# plot(st_geometry(world %>% filter(continent == "Africa")), add=TRUE)

writeRaster(ras_mask,
            "data/stable_transmission_mask.grd",
            overwrite = TRUE)
