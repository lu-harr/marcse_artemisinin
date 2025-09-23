# coverage probabilities: what's the probability a prediction interval covers 
# the true parameter? (Which I'm interpreting to mean the probability that the 
# prediction interval covers the realised data - check with JF/NG)
# or perhaps I want the probability that data (with associated uncertainty) covers 
# expectation (model estimate)?

setwd("~/Desktop/MARCSE/k13_seafrica")
source("code/setup.R")
source("code/vis/vis_funcs.R")
library(sf)

preds <- rast(paste0("output/", "k13_marcse", "/gneiting_ahmc/preds_all.tif")) %>%
  aggregate(fact=3)

# for all points in prevalence dataset,
# does the CI cover the prevalence at the point?
# (imagine we struggle at zero?)

rast_plot(preds[[which(grepl("CI", names(preds)))]])
ggsave("~/Desktop/ci_width.png")

preds <- rast("output/k13_marcse/gneiting_ahmc/preds_upper_lower.tif")
names(preds)
dat <- setup_mut_data("data/clean/moldm_marcse_k13_nomarker.csv")

# "are the data inside the confidence intervals ?"
coverage_prob <- function(preds, dat){
  dat <- dat %>%
    mutate(prevalence = present / tested,
           upper = NA,
           lower = NA)
  yrs <- unique(dat$year)
  for (yr in yrs){
    idx <- which(dat$year == yr)
    tmp <- subset(preds, grepl(yr, names(preds)))
    lower <- terra::extract(subset(tmp, grepl("2.5", names(tmp))), dat[idx, c("x", "y")], ID=FALSE)
    upper <- terra::extract(subset(tmp, grepl("97.5", names(tmp))), dat[idx, c("x", "y")], ID=FALSE)
    dat[idx, c("upper", "lower")] <- cbind(upper, lower)
  }
  
  dat$covered <- dat$prevalence >= dat$lower & dat$prevalence <= dat$upper
  sum(dat$covered, na.rm=TRUE) / nrow(dat)
  sum(dat$covered[dat$prevalence > 0], na.rm=TRUE) / sum(dat$prevalence > 0, na.rm=TRUE)
  
  dat %>% filter(prevalence > 0)
}







