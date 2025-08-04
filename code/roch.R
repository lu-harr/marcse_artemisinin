preds <- rast("output/k13_marcse/gneiting_sparse/preds_all.grd") 
ras <- rast("output/k13_marcse/surveillance_effort_k13_marcse.grd")

# does uncertainty look anything like rate of change?
plot(preds$`2022_post_sd`)
plot(preds$`2022_post_median`)
#focalMat(preds, 0.05, "Gauss")

sample_var <- function(x){sum((x - mean(x, na.rm = TRUE))**2, na.rm = TRUE) / (sum(!is.na(x)))}

tmp = focal(preds$`2022_post_median`, 7, sample_var)

plot(sqrt(tmp), main = "Local sample std dev")


preds <- rast("output/mdr1246/gneiting_sparse/preds_all.grd")

ras <- rast("output/mdr1246/surveillance_effort_mdr1246_annual.grd") %>%
  tidyterra::select(surv_2006, surv_2005, surv_2007) %>%
  mask(preds$`2000_post_median`) %>%
  trim() %>%
  sum()

preds <- trim(preds)

tmp = focal(preds$`2006_post_median`, 7, sample_var) %>%
  mask(preds$`2000_post_median`)

{png("~/Desktop/presentations/marcse/uncertainty.png", 
     height=1600, width=2200, pointsize=40)
par(mfrow = c(2, 3), oma = c(0,0,2,0))
plot(preds$`2006_post_median`, main = "Estimated prevalence")
plot(preds$`2006_post_sd`, main = "Posterior std deviation")
plot(preds$`2006_post_median`, preds$`2006_post_sd`, xlab = "Posterior median", ylab = "Posterior SD")
plot(sqrt(tmp), main = "Local sample std dev")
plot(ras$sum, main = "Surveillance effort")
mtext("Uncertainty case study: Pfmdr1 D1246Y, 2006", outer = TRUE, font = 2)
dev.off()}


