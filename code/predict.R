source("code/setup.R")

mut_data <- read_rds("output/mat52_model/mut_data.rds")
parameters <- read_rds("output/mat52_model/parameters.rds")
# read_rds(mut_prob_obs, "output/mat52_model/mut_prob_obs.rds")
kernel <- read_rds("output/mat52_model/kernel.rds")
random_field <- read_rds("output/mat52_model/random_field.rds")
m <- read_rds("output/mat52_model/m.rds")
draws <- read_rds("output/mat52_model/draws.rds")

# mut_data <- read_rds("../dhps_mapping/quintuple_question/quintuple/mut_data.rds")
# parameters <- read_rds("../dhps_mapping/quintuple_question/quintuple/parameters.rds")
# # read_rds(mut_prob_obs, "../dhps_mapping/quintuple_question/quintuple/mut_prob_obs.rds")
# kernel <- read_rds("../dhps_mapping/quintuple_question/quintuple/kernel.rds")
# random_field <- read_rds("../dhps_mapping/quintuple_question/quintuple/random_field.rds")
# m <- read_rds("../dhps_mapping/quintuple_question/quintuple/m.rds")
# draws <- read_rds("../dhps_mapping/quintuple_question/quintuple/draws.rds")
# 
# mut_data <- read_rds("../dhps_mapping/quintuple_question/540E/mut_data.rds")
# parameters <- read_rds("../dhps_mapping/quintuple_question/540E/parameters.rds")
# # read_rds(mut_prob_obs, "../dhps_mapping/quintuple_question/540E/mut_prob_obs.rds")
# kernel <- read_rds("../dhps_mapping/quintuple_question/540E/kernel.rds")
# random_field <- read_rds("../dhps_mapping/quintuple_question/540E/random_field.rds")
# m <- read_rds("../dhps_mapping/quintuple_question/540E/m.rds")
# draws <- read_rds("../dhps_mapping/quintuple_question/540E/draws.rds")

source("code/build_design_matrix.R")

# dropping prediction code from dhps_africa_greta
stable_transmission_mask <- covariates$pfpr_2010 %>%
  aggregate(15)
stable_transmission_mask[stable_transmission_mask < -0.62] <- NA

source("code/predict_to_raster.R")

year1 = 2000
year2 = 2010
year3 = 2020


tmp2 <- predict_to_ras(covariates,
                       year1,
                       draws,
                       parameters,
                       agg_factor = 15,
                       stable_transmission_mask = stable_transmission_mask)

tmp3 <- predict_to_ras(covariates,
                       year2,
                       draws,
                       parameters,
                       agg_factor = 15,
                       stable_transmission_mask = stable_transmission_mask)

tmp4 <- predict_to_ras(covariates,
                       year3,
                       draws,
                       parameters,
                       agg_factor = 15,
                       stable_transmission_mask = stable_transmission_mask)

library(raster)
tmp2 <- stack(tmp2)
tmp3 <- stack(tmp3)
tmp4 <- stack(tmp4)

writeRaster(tmp2, paste0("output/mat52_model/k13_", year1, ".tif"), overwrite=TRUE)
writeRaster(tmp3, paste0("output/mat52_model/k13_", year2, ".tif"), overwrite=TRUE)
writeRaster(tmp4, paste0("output/mat52_model/k13_", year3, ".tif"), overwrite=TRUE)

library(viridisLite)

{png("figures/k13.png",
     height=2800, width=2800, pointsize=40)
  par(mfrow=c(3,3), mar=c(0.2, 0, 0, 5.5), oma=c(4,4.5,3,0))
  plot(sqrt(trim(tmp2$post_mean)), col="grey90", legend=F, xaxt="n")
  tmp_pts <- mut_data %>%
    filter(year == 2009)
  points(tmp_pts[,c("x", "y")],
         col=ifelse(tmp_pts$present == 0, "grey50", "orange"), lwd=3, cex=1.2)
  mtext(year1, line=3, side=2)
  mtext("Data", line=1)
  
  plot(sqrt(trim(tmp2$post_mean)), col=viridis(100),
       breaks=seq(0, 0.6, length.out=100), legend=F, xaxt="n", yaxt="n")
  mtext("Mean", line=1)
  
  plot(sqrt(trim(tmp2$post_sd)), col=viridis(100), xaxt="n", yaxt="n",
       breaks=seq(0, 0.34, length.out=100), legend=F)
  mtext("SD", line=1)
  
  
  plot(sqrt(trim(tmp2$post_mean)), col="grey90", legend=F)
  tmp_pts <- mut_data %>%
    filter(year == 2019)
  points(tmp_pts[,c("x", "y")],
         col = ifelse(tmp_pts$present == 0, "grey50", "orange"), lwd=2,
         cex = 1.2 + tmp_pts$present/tmp_pts$tested * 10)
  mtext(year2, line=3, side=2)
  
  plot(sqrt(trim(tmp3$post_mean)), col=viridis(100),
       breaks=seq(0,0.6, length.out=100), legend=F, yaxt="n")
  
  plot(sqrt(trim(tmp3$post_sd)), col=viridis(100), yaxt="n",
       breaks=seq(0, 0.34, length.out=100), legend=F)
  #mtext("k13 prevalence: preliminary model", outer=TRUE)
  
  plot(sqrt(trim(tmp4$post_mean)), col="grey90", legend=F)
  tmp_pts <- mut_data %>%
    filter(year > 2021)
  points(tmp_pts[,c("x", "y")],
         col = ifelse(tmp_pts$present == 0, "grey50", "orange"), lwd=2,
         cex = 1.2 + tmp_pts$present/tmp_pts$tested * 10)
  mtext(year3, line=3, side=2)
  
  plot(sqrt(trim(tmp4$post_mean)), col=viridis(100),
       breaks=seq(0,0.6, length.out=100), legend=F, yaxt="n")
  
  plot(sqrt(trim(tmp4$post_sd)), col=viridis(100), yaxt="n",
       breaks=seq(0, 0.34, length.out=100), legend=F)
  
  legend_tix <- c(0, 0.1, 0.2, 0.3)
  par(new=TRUE, mfrow=c(1,3), mfg=c(1,2))
  plot(0, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
  plot(sqrt(trim(tmp2$post_mean)), col=viridis(100),
       breaks=seq(0, 0.6, length.out=100), xaxt="n", yaxt="n", legend.only=TRUE,
       axis.args=list(at = sqrt(legend_tix), labels = legend_tix), 
       legend.width=1.2)
  
  legend_tix <- c(0, 0.025, 0.05, 0.075, 0.1)
  plot(0, type="n", xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
  legend("center", c("Presence", "Absence"), fill=c("orange", "grey50"))
  plot(sqrt(trim(tmp2$post_sd)), col=viridis(100),
       breaks=seq(0, 0.34, length.out=100), xaxt="n", yaxt="n", legend.only=TRUE,
       axis.args=list(at = sqrt(legend_tix), labels = legend_tix),
       legend.width=1.2)
  
  
  dev.off()}

all_outs <- stack(tmp2, tmp3, tmp4) %>%
  getValues() %>%
  as.data.frame() %>%
  cbind(xyFromCell(tmp2, seq_len(ncell(tmp2)))) %>%
  setNames(c("Mean 2010", "SD 2010", "Mean 2017", "SD 2017", 
             "Mean 2022", "SD 2022", "x", "y"))

st = stack(tmp2, tmp3, tmp4)
names(st) <- c("Mean 2010", "SD 2010", "Mean 2017", "SD 2017", 
               "Mean 2022", "SD 2022")
coords <- xyFromCell(st, seq_len(ncell(st)))
st <- stack(as.data.frame(getValues(st)))
names(st) <- c("value", "variable")
st <- mutate(st, variable = gsub("\\.", " ", variable))
  
st <- cbind(coords, st)

st_mean <- st %>%
  filter(grepl("Mean", variable)) %>%
  mutate(variable = gsub("Mean ", "", variable))

st_sd <- st %>%
  filter(grepl("SD", variable)) %>%
  mutate(variable = gsub("SD ", "", variable))

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(dplyr)
library(terra)
world <- ne_countries(scale="medium", returnclass = "sf")

ggplot() +
  geom_sf(data = filter(world, continent == "Africa"), fill = "white") +
  geom_tile(data = st_mean, aes(x, y, fill = value)) +
  facet_wrap(~ factor(variable, c("2010", "2017", "2022")), nrow=1) +
  scale_fill_viridis_c(na.value = NA, "Prevalence") +
  xlab("Longitude") +
  ylab("Latitude")
ggsave("k13_mean.png", width = 6, height = 2, scale = 2)

ggplot() +
  geom_sf(data = filter(world, continent == "Africa"), fill = "white") +
  geom_tile(data = st_sd, aes(x, y, fill = value)) +
  facet_wrap(~ factor(variable, c("2010", "2017", "2022")), nrow=1) +
  scale_fill_distiller(palette="YlOrBr", na.value = NA, "Uncertainty", trans="sqrt") +
  xlab("Longitude") +
  ylab("Latitude")
ggsave("k13_sd.png", width = 6, height = 2, scale = 2)

# perhaps what I should do is give up on the idea that they should be one plot


# library(rasterVis)
# gplot(stack(tmp2, tmp3, tmp4)) +
#   geom_tile(aes(fill = value)) +
#   facet_wrap(~ variable, ncol = 2) +
#   scale_fill_viridis_c(na.value = "white") +
#   coord_equal()








