# cook a stable transmission mask ... 

source("code/setup.R")
library(sf)

ras <- covariates$pfpr_2022 #%>%
  #aggregate(AGG_FACTOR)

# was initially applying some blurring to get a smoother surface:
# wt <- focalMat(ras, 0.5, type = "Gauss")
# ras_blur <- focal(ras - minmax(ras)[1], wt, na.rm=TRUE)

# ras_mask <- ras
# ras_mask[ras_mask < 0.01] <- NA
# ras_mask <- mask(ras_mask, ras)
# ras_mask <- mask(ras_mask, 
#                  world %>% 
#                    filter(continent == "Africa") %>% 
#                    st_buffer(dist = 100000))

# I'm not sure why this is negative ... on MAP's dashboard they've scaled it from 
# a "rate" (they're definitely still calling it a rate) to [0,100] (% prevalence 
# in 2-10 yos)

# This raster looks a bit more conservative in the south, e.g. Namibia, Botswana:
# k13 has popped up there but that info isn't in the Surveyor (right?)

ras_mask <- ras
# this is about as much as we're going to squeeze out of this ...
ras_mask[ras_mask < -0.726] <- NA

plot(ras_mask)

library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale="medium", returnclass = "sf")

# country shps for plotting/masking
ssafr <- world %>%
  filter(continent == "Africa") %>%
  filter(!name %in% c("Algeria", "Libya", "Egypt", "Tunisia", "Morocco")) %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37))

ras_mask <- mask(ras_mask,
                 ssafr)
plot(ras_mask)
plot(st_geometry(st_as_sf(ssafr)), add=TRUE)

writeRaster(ras_mask,
            "data/stable_transmission_mask.grd",
            overwrite = TRUE)


ras <- rast("data/stable_transmission_mask.grd")
plot(ras)

