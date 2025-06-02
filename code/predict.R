mut_data <- read_rds("output/mat52_model/mut_data.rds")
parameters <- read_rds("output/mat52_model/parameters.rds")
# read_rds(mut_prob_obs, "output/mat52_model/mut_prob_obs.rds")
kernel <- read_rds("output/mat52_model/kernel.rds")
random_field <- read_rds("output/mat52_model/random_field.rds")
m <- read_rds("output/mat52_model/m.rds")
draws <- read_rds("output/mat52_model/draws.rds")

source("code/setup.R")
source("code/build_design_matrix.R")
source("code/predict_to_raster.R")

# dropping prediction code from dhps_africa_greta
stable_transmission_mask <- covariates$pfpr_2010 %>%
  aggregate(15)
stable_transmission_mask[stable_transmission_mask < -0.62] <- NA
# stable_transmission_mask <- crop(stable_transmission_mask,
#                                  ext(-21, 60, -35.0833320617676, 37.4166679382324))

tmp2 <- predict_to_ras(covariates,
                       2010,
                       draws,
                       parameters,
                       agg_factor = 15,
                       stable_transmission_mask = stable_transmission_mask)

tmp3 <- predict_to_ras(covariates,
                       2020,
                       draws,
                       parameters,
                       agg_factor = 15,
                       stable_transmission_mask = stable_transmission_mask)

tmp4 <- predict_to_ras(covariates,
                       2022,
                       draws,
                       parameters,
                       agg_factor = 15,
                       stable_transmission_mask = stable_transmission_mask)

library(raster)
tmp2 <- stack(tmp2)
tmp3 <- stack(tmp3)
tmp4 <- stack(tmp4)

{png("figures/k13.png",
     height=2800, width=2800, pointsize=40)
  par(mfrow=c(3,3), mar=c(0.2, 0, 0, 5.5), oma=c(4,4.5,3,0))
  plot(sqrt(trim(tmp2$post_mean)), col="grey90", legend=F, xaxt="n")
  tmp_pts <- mut_data %>%
    filter(year == 2009)
  points(tmp_pts[,c("x", "y")],
         col=ifelse(tmp_pts$present == 0, "grey50", "orange"), lwd=3, cex=1.2)
  mtext("2009", line=3, side=2)
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
  mtext("2019", line=3, side=2)
  
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
  mtext("2022", line=3, side=2)
  
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
