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

preds <- rast("output/mdr184/gneiting_sparse/preds_all.grd") 
# does uncertainty look anything like rate of change?
plot(preds$`2022_post_sd`)
plot(preds$`2022_post_median`)
plot(preds$`2022_post_median`, preds$`2022_post_sd`)

df <- gg_ras_prep(preds)$df

df <- pivot_wider(df, names_from = tag, values_from = val, id_cols = c(x, y, year))

tmp2 <- slice_sample(df, n = 1e5)

tmp2 = df %>% filter(year == 2022)

ggplot(tmp2) +
  geom_point(aes(x = medi, y = sd),
            col = "darkgrey", alpha = 0.2, size = 0.8) +
  xlab("Median") +
  ylab("Std Deviation") +
  theme_bw() #+
  #facet_wrap(~year)
ggsave("~/Desktop/presentations/MARCSE/mdr184_sd_over_medi.png", 
       height = 5, width = 5, scale = 0.5)

preds <- rast("output/mdr86/gneiting_sparse/preds_all.grd") 
preds <- rast(list(preds$`2006_post_median`,
              preds$`2006_post_sd`))
# idk why this now ain't working
df <- gg_ras_prep(preds)$df %>%
  pivot_wider(names_from = tag, values_from = val, id_cols = c(x, y, year))

ggplot(df) +
  geom_point(aes(x = medi, y = sd),
             col = "darkgrey", alpha = 0.2, size = 0.8) +
  xlab("Median") +
  ylab("Std Deviation") +
  theme_bw()
ggsave("~/Desktop/presentations/MARCSE/mdr86_sd_over_medi.png", 
       height = 5, width = 5, scale = 0.5)

preds <- rast("output/crt76/gneiting_sparse/preds_all.grd") 
preds <- rast(list(preds$`2014_post_median`,
                   preds$`2014_post_sd`))
df <- gg_ras_prep(preds)$df %>%
  pivot_wider(names_from = tag, values_from = val, id_cols = c(x, y, year))

ggplot(df) +
  geom_point(aes(x = medi, y = sd),
             col = "darkgrey", alpha = 0.2, size = 0.8) +
  xlab("Median") +
  ylab("Std Deviation") +
  theme_bw()
ggsave("~/Desktop/presentations/MARCSE/crt76_sd_over_medi.png", 
       height = 5, width = 5, scale = 0.5)



