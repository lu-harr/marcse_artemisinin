# folding functions for plotting in here
library(looseVis)
library(cowplot)
library(iddoPal)
library(sf)

obs_prev_panel <- function(data_path, 
                           pred_path, 
                           main = "", 
                           show_nas = FALSE, 
                           #pal = colorRamp(viridis(10)), 
                           pal = colorRamp(iddo_palettes$BlGyRd),
                           xlim = c(0,1), # define limits to pred/obs panel
                           ylim = c(0,1), # define limits to pred/obs panel
                           facet_bins = NULL, # apply facets over time?
                           ave_tag = "_50", # mean? median? what are the surfaces called in the stack?
                           buffer = 1){
  # 
  mut_data <- setup_mut_data(data_path, min_year = 2000)
  preds <- rast(pred_path)
  
  # get predictions for each row in `mut_data`
  mut_data$pred <- NA
  yrs_to_extract <- unique(mut_data$year)
  for (yr in yrs_to_extract){
    if (yr %in% pfpr_years){
      idx <- which(mut_data$year == yr)
      val <- terra::extract(preds[[paste0(yr, ave_tag)]], 
                            mut_data[idx, c("x", "y")],
                            ID = FALSE, search_radius = buffer)
      if(ncol(val) < 3){
        # idk why we have to have an inconsistent return when |idx| == 1
        mut_data[idx, "pred"] <- val[1,1]}
      else{
        mut_data[idx, "pred"] <- val[, paste0(yr, ave_tag)]
      }
    }
  }
  
  mut_data$cex <- scale_cex(mut_data$tested, sqrt, max_cex = 5)
  un_pred <- mut_data[is.na(mut_data$pred),]
  mut_data <- mut_data[!is.na(mut_data$pred),]
  
  # mut_data$diffs = abs(mut_data$present / mut_data$tested - mut_data$pred)
  mut_data$diffs <- mut_data$present / mut_data$tested - mut_data$pred
  mut_data <- arrange(mut_data, abs(diffs))
  mut_data$diffs_scaled = mut_data$diffs / 2 + 0.5 # hopefully grey ends up where diffs == 0?
  
  if (!is.null(facet_bins)){
    mut_data$year_bin <- cut(mut_data$year, 
                             c(min(mut_data$year) - 1, facet_bins, max(mut_data$year)))
  }
  
  p1 <- ggplot(mut_data) +
    geom_point(mapping = aes(x = present / tested, y = pred), 
               size = mut_data$cex,
               col = rgb(pal(mut_data$diffs_scaled), maxColorValue = 255),
               shape = 1) +
    geom_abline(slope = 1, intercept = 0) +
    xlim(xlim) +
    ylim(ylim) +
    xlab("Observed prevalence") +
    ylab("Predicted prevalence") +
    theme_bw() #+
    #theme(plot.margin = unit(c(0.2,0.2,2.1,0.2), "cm"))
  
  p2 <- ggplot(mut_data) +
    geom_point(mapping = aes(x = year, y = present / tested - pred), 
               size = mut_data$cex,
               col = rgb(pal(mut_data$diffs_scaled), maxColorValue = 255),
               shape = 1) +
    geom_hline(yintercept = 0) +
    xlab("Observed prevalence") +
    ylab("Predicted prevalence") +
    theme_bw()
  
  # this is a bit hacky but I want to constrain the endpoints of my colour scale
  # so that they mean roughly the same between different markers
  cols = seq(min(mut_data$diffs_scaled), 
             max(mut_data$diffs_scaled), length.out = 100) %>%
    pal() %>%
    rgb(maxColorValue = 255)
  
  p3 <- ggplot() +
    geom_sf(data = st_as_sf(afr), fill = "white") + # not showing anything in the background here ...
    geom_point(data = mut_data, 
               aes(x = x, y = y, col = diffs),
               shape = 1, size = mut_data$cex) +
    scale_color_gradientn(colours = cols, 
                          name = "Residuals\n(Observed - Predicted)") + # this needs re-scaling (back to what it was..)
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    labs(title = main) +
    theme(legend.position = "bottom")
  
  if (show_nas){ # not sure how this plays with faceting
    p3 <- p3 + geom_point(data = un_pred, 
                          aes(x = x, y = y), 
                          col = "orange", size = un_pred$cex, shape = 1)
  }
  
  if (!is.null(facet_bins)){
    p1 <- p1 + facet_wrap(vars(.data$year_bin), ncol = 1)
    p3 <- p3 + facet_wrap(vars(.data$year_bin), ncol = 1) 
  }
  
  #plot_grid(p1, p3, rel_widths = c(0.5, 0.52))
  
  plot_grid(p1, p2, ncol = 1) %>%
    plot_grid(p3, rel_widths = c(0.4, 0.6))
}



# e.g.:
# might want to re-land some points inside of model fitting
obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
               "output/k13_marcse/gneiting_sparse/preds_all.grd",
               #main = "k13 gneiting binom", 
               xlim = c(0, 0.6), ylim = c(0, 0.6),
               ave_tag = "_post_median", buffer = 100000)
ggsave("~/Desktop/presentations/marcse/residuals_k13m_t.png", height = 3.7, width = 5, scale = 1.5)

# could try giving it longer ?
obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
               "output/k13_marcse/bb_gne/preds_all.tif",
               main = "k13 betabinom gneiting", 
               xlim = c(0, 0.6), ylim = c(0, 0.6),
               show_nas = TRUE, buffer = 100000)
# bit spooked by the points changing between these two ...
# might be points falling off the mask?
# that is so many points !








