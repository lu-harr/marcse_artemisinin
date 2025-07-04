# cook a stable transmission mask ... 
# subject to change but would like to be smoothish
# (get feedback from Jen/Philippe/Karen before publication obvs)
source("code/setup.R")
library(sf)

ras <- covariates$pfpr_2022 #%>%
  #aggregate(AGG_FACTOR)

wt <- focalMat(ras, 0.5, type = "Gauss")
ras_blur <- focal(ras - minmax(ras)[1], wt, na.rm=TRUE)

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
ras_mask[ras_mask < -0.65] <- NA

plot(ras_mask)

library(rnaturalearth)
library(rnaturalearthdata)
world <- ne_countries(scale="medium", returnclass = "sf")

# country shps for plotting/masking
afr <- world %>%
  filter(continent == "Africa") %>%
  filter(!name %in% c("Algeria", "Libya", "Egypt", "Tunisia", "Morocco")) %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37))

ras_mask <- mask(ras_mask,
                 afr)

writeRaster(ras_mask,
            "data/stable_transmission_mask.grd",
            overwrite = TRUE)





