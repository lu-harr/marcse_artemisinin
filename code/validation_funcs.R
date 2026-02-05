# fucntions for validation scripts

extract_preds <- function(data_path,
                          pred_path,
                          ave_tag = "_50",
                          buffer = 0){
  # grab predictions for all points in dataset for un-held-out models
  # bring in coords associated with predictions
  mut_data <- setup_mut_data(data_path, min_year = MIN_YEAR, buffer = buffer)
  preds <- rast(pred_path)
  yrs_pred <- str_extract(names(preds), "\\d{4}")
  
  # get predictions for each row in `mut_data`
  mut_data$pred <- NA
  yrs_to_extract <- unique(mut_data$year)
  for (yr in yrs_to_extract){
    message(yr)
    if (yr %in% yrs_pred){
      idx <- which(mut_data$year == yr)
      val <- terra::extract(preds[[paste0(yr, ave_tag)]], 
                            mut_data[idx, c("x", "y")],
                            ID = FALSE, 
                            search_radius = buffer)
      
      if (length(idx) == 1){
        # not quite sure why, but we get a 2*2 with predicted value twice ...
        val <- val[1, ]
      }
      
      mut_data[idx, "pred"] <- val[, paste0(yr, ave_tag)]
    #plot(mut_data$present[idx] / mut_data$tested[idx], mut_data$pred[idx], main = yr)
    }
  }
  
  mut_data
}

# e.g.:
# mut_data <- extract_preds(data_path = "data/clean/moldm_marcse_k13_nomarker.csv",
#                           pred_path = "output/k13_marcse/bb_gne/preds_medians.tif")
# plot(mut_data$present / mut_data$tested, mut_data$pred)

extract_preds_cv <- function(data_path,
                              pred_path,
                              folds,
                              ave_tag = "_50", # vestigial at this point
                              in_buffer = 0,
                              out_buffer = 0 # integer for dragging points back to land
                              ){
  # grab predictions for all points in dataset on held-out models

  folds <- lapply(1:length(folds), function(idx){
    data.frame(fold = idx,
               pts = folds[[idx]])
  }) %>%
    do.call(bind_rows, .) %>%
    arrange(pts)
  
  message(paste(unique(folds$fold)))
  
  # bring in coords associated with predictions + associate with fold
  mut_data <- setup_mut_data(data_path, min_year = MIN_YEAR, buffer = in_buffer) %>%
    mutate(fold = folds$fold)
  
  message(nrow(mut_data))
  
  
  preds_avail <- list.files(pred_path)
  preds_avail <- preds_avail[grep("preds_medians", preds_avail)]
  # (in case not all folds predicted to:)
  folds_avail <- unique(str_extract(preds_avail, "(?<=_)\\d+(?=\\.tif)"))
  message(preds_avail)
  message(folds_avail) # am I missing one?
  
  # get predictions for each row in `mut_data`
  mut_data$pred <- NA
  extract_attempted = c()
  
  for (fidx in 1:length(folds_avail)){
    preds <- rast(paste0(pred_path, preds_avail[fidx]))
    yrs_pred <- str_extract(names(preds), "\\d{4}")
    
    message(paste0(fidx, "!", folds_avail[fidx])) # something happened to 2: no yrs_to_extract
    
    yrs_to_extract <- mut_data %>%
      filter(fold == folds_avail[fidx]) %>%
      dplyr::select(year) %>%
      unique() %>%
      unlist()
    
    message(yrs_to_extract)
    
    for (yr in yrs_to_extract){
      if (yr %in% yrs_pred){
        idx <- which(mut_data$year == yr & mut_data$fold == folds_avail[fidx])
        extract_attempted  <- c(extract_attempted, idx)
        
        val <- terra::extract(preds[[paste0(yr, ave_tag)]], 
                              mut_data[idx, c("x", "y")],
                              ID = FALSE, search_radius = out_buffer)
        
        message(paste(length(idx), ",", sum(is.na(val))))
  
        mut_data[idx, "pred"] <- val[, paste0(yr, ave_tag)]
      } else {
        message(paste(yr, " not in preds"))
      }
    }
  }
  # is the problem that I'm not extracting for some points or that extraction
  # is off-land?
  message(length(setdiff(1:nrow(mut_data), extract_attempted)))


  mut_data
}

# folds <- read_rds("output/k13_marcse/bb_gne/cv_folds.rds")

# tmp = extract_preds_cv(data_path = "data/clean/moldm_marcse_k13_nomarker.csv",
#                  pred_path = "output/k13_marcse/bb_gne/cv_preds/",
#                  folds = folds)


rmse <- function(mut_data){
  # root mean square error on probability scale
  sqrt(sum((mut_data$pred - mut_data$present/mut_data$tested)**2, na.rm = TRUE) / 
         sum(!is.na(mut_data$pred)))
}

unadjusted_rsq <- function(dat){
  tmp <- filter(dat, !is.na(pred))
  preds <- tmp$pred
  actual <- tmp$present/tmp$tested
  rss <- sum((preds - actual) ^ 2)
  tss <- sum((actual - mean(actual)) ^ 2)
  1 - rss/tss
}

# or using the built-in:
# unadjusted_rsq <- function (x, y) cor(x, y) ^ 2

nn_measure <- function(mut_data, draws_path){
  # from YSF: include nearest neighbour measure/ some other proximity measure
  # in validation panels to understand if mis-prediction is caused by multiple
  # observations near each other with variable measurements
  # need to have years scaled
  scaled_years <- scale_years(range(mut_data$year))
  mut_data$year_scaled <- scaled_years[as.character(mut_data$year)]
  
  mut_data$xr <- degrees_to_radians(mut_data$x)
  mut_data$yr <- degrees_to_radians(mut_data$y)
  
  # get dims right for tf_gneiting
  mut_data_nume <- mut_data %>%
    mutate(year_scaled = as.numeric(year_scaled)) %>%
    dplyr::select(xr, yr, year_scaled) %>%
    as.matrix() %>%
    tf$expand_dims(0L)
  
  draws <- readRDS(paste0(draws_path, "draws.rds")) %>%
    summary()
  params <- draws$quantiles[c("gneiting_len", "gneiting_tim", "gneiting_sd"), "50%"]
  
  covs <- tf_gneiting(mut_data_nume,
              mut_data_nume,
              lengthscale = params[1], 
              timescale = params[2],
              variance = params[3]**2, 
              active_dims = c(0,1,2)) %>%
    tf$squeeze() %>%
    as.matrix()
  
  # YSF has the maximum of the covs for a point compared to *all* training points
  # yet to implement folds, so for now, find max for all points
  diag(covs) = NA # hack - diag won't be in training set
  apply(covs, 1, mean, na.rm=TRUE)
}



obs_prev_panel <- function(mut_dat, 
                           main = "", 
                           show_nas = FALSE, 
                           pal = colorRamp(iddo_palettes$BlGyRd),
                           xlim = c(0,1), # define limits to pred/obs panel
                           ylim = c(0,1), # define limits to pred/obs panel
                           facet_bins = NULL, # apply facets over time?
                           ave_tag = "_50", # mean? median? what are the surfaces called in the stack?
                           buffer = 1, # option to reland points?
                           bb = NULL, # bounding box for inset
                           as_row = FALSE){
  # observed vs predicted values: give me a three-panelled plot
  
  mut_dat$cex <- scale_cex(mut_dat$tested, sqrt, max_cex = 5)
  un_pred <- mut_dat[is.na(mut_dat$pred),]
  mut_dat <- mut_dat[!is.na(mut_dat$pred),]
  
  mut_dat$diffs <- mut_dat$present / mut_dat$tested - mut_dat$pred
  mut_dat <- arrange(mut_dat, abs(diffs))
  
  # grey should end up where diffs == 0:
  diffs_ext <- max(abs(mut_dat$diffs), na.rm=TRUE)
  mut_dat$diffs_scaled = mut_dat$diffs / 2 + 0.5
  mut_dat$diffs_scaled = mut_dat$diffs / (diffs_ext * 2) + 0.5
  
  if (!is.null(facet_bins)){
    mut_dat$year_bin <- cut(mut_dat$year, 
                             c(min(mut_dat$year) - 1, facet_bins, max(mut_dat$year)))
  }
  
  p1 <- ggplot(mut_dat) +
    geom_point(mapping = aes(x = present / tested, y = pred), 
               size = mut_dat$cex,
               col = rgb(pal(mut_dat$diffs_scaled), maxColorValue = 255),
               shape = 1) +
    geom_abline(slope = 1, intercept = 0) +
    xlim(xlim) +
    ylim(ylim) +
    xlab("Observed prevalence") +
    ylab("Predicted prevalence") +
    theme_bw() #+
  #theme(plot.margin = unit(c(0.2,0.2,2.1,0.2), "cm"))
  
  p2 <- ggplot(mut_dat) +
    geom_point(mapping = aes(x = year, y = present / tested - pred), 
               size = mut_dat$cex,
               col = rgb(pal(mut_dat$diffs_scaled), maxColorValue = 255),
               shape = 1) +
    geom_hline(yintercept = 0) +
    xlab("Year") +
    ylab("Observed - Predicted prevalence") +
    theme_bw()
  
  # this is a bit hacky but I want to constrain the endpoints of my colour scale
  # so that they mean roughly the same between different markers
  cols = seq(min(mut_dat$diffs_scaled, na.rm=T), 
             max(mut_dat$diffs_scaled, na.rm=T), length.out = 100) %>%
    pal() %>%
    rgb(maxColorValue = 255)
  
  p3 <- ggplot() +
    geom_sf(data = st_as_sf(afr), fill = "white") + # not showing anything in the background here ...
    geom_point(data = mut_dat, 
               aes(x = x, y = y, col = diffs),
               shape = 1, size = mut_dat$cex) +
    scale_color_gradientn(colours = cols, 
                          name = "Residuals\n(Observed - Predicted)") + 
    # this needs re-scaling (back to what it was..)
    theme_bw() +
    xlab("Longitude") +
    ylab("Latitude") +
    theme(legend.position = "bottom")
  
  if (!is.null(bb)){
    # inset panel: RWA/UGA
    p4 <- p3 +
      xlim(bb[1:2]) +
      ylim(bb[3:4]) +
      theme(axis.title = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.position = "none",
            plot.margin = margin(0, 0, 0, 0)) %>%
      suppressWarnings()
    
    p3 <- ggdraw() +
      draw_plot(p3 +
                  geom_rect(aes(xmin = bb[1], xmax = bb[2], 
                                ymin = bb[3], ymax = bb[4]),
                            linetype = "dashed",
                            fill = NA, col = "grey30", lwd=0.2)) +
      draw_plot(p4, x = 0.15, y = 0.21, width = 0.3, height = 0.3)%>%
      suppressWarnings()
    #draw_plot(p4, x = -20, y = -30, width = 50, height = 50)
  }
  
  if (show_nas){ # not sure how this plays with faceting
    p3 <- p3 + geom_point(data = un_pred, 
                          aes(x = x, y = y), 
                          col = "orange", size = un_pred$cex, shape = 1)
  }
  
  if (!is.null(facet_bins)){
    p1 <- p1 + facet_wrap(vars(.data$year_bin), ncol = 1)
    p3 <- p3 + facet_wrap(vars(.data$year_bin), ncol = 1) 
  }
  
  if (!as_row){
    plot_grid(p1, p2, ncol = 1) %>%
      plot_grid(p3, rel_widths = c(0.4, 0.7))
  } else {
    plot_grid(p1, p2, p3 + theme(legend.position = "none"), nrow = 1) %>%
      plot_grid(get_plot_component(p3, "guide-box", return_all = TRUE)[[3]], 
                ncol = 1, rel_heights = c(1, 0.3))
  }
}



obs_prev_panel_nn <- function(data_path, 
                           pred_path, 
                           draws_path,
                           main = "", 
                           show_nas = FALSE, 
                           pal = colorRamp(iddo_palettes$BlGyRd),
                           xlim = c(0,1), # define limits to pred/obs panel
                           ylim = c(0,1), # define limits to pred/obs panel
                           facet_bins = NULL, # apply facets over time?
                           ave_tag = "_50", # mean? median? what are the surfaces called in the stack?
                           buffer = 1 # option to reland points?
                           ){
  # visualise with NN metric
  
  mut_data <- extract_preds(data_path, pred_path, ave_tag, buffer)
  mut_data$nn <- nn_measure(mut_data, draws_path)
  mut_data$cex <- scale_cex(mut_data$tested, sqrt, max_cex = 5)
  un_pred <- mut_data[is.na(mut_data$pred),]
  mut_data <- mut_data[!is.na(mut_data$pred),]
  
  message(paste0("Sites with preds: ", nrow(mut_data)))
  message(paste0("Mean error: ", mean(mut_data$pred - mut_data$present/mut_data$tested)))
  message(paste0("Mean abs error: ", mean(abs(mut_data$pred - mut_data$present/mut_data$tested))))
  
  mut_data$diffs <- mut_data$present / mut_data$tested - mut_data$pred
  mut_data <- arrange(mut_data, abs(diffs))
  
  p1 <- ggplot(mut_data) +
    geom_point(mapping = aes(x = present / tested, y = pred, col = nn), 
               size = mut_data$cex,
               shape = 1) +
    geom_abline(slope = 1, intercept = 0) +
    xlim(xlim) +
    ylim(ylim) +
    xlab("Observed prevalence") +
    ylab("Predicted prevalence") +
    theme_bw()
  
  p2 <- ggplot(mut_data) +
    geom_point(mapping = aes(x = year, y = present / tested - pred, col = nn), 
               size = mut_data$cex,
               shape = 1) +
    geom_hline(yintercept = 0) +
    xlab("Year") +
    ylab("Observed - Predicted prevalence") +
    theme_bw()

  plot_grid(p1, p2, nrow = 1)
}




coverage_probabilities_from_observation_model <- function(mut_data,
                                                         draws_path,
                                                         nsim = 100,
                                                         probs){
  # given draws, grab median rho
  draws <- readRDS(paste0(draws_path, "draws.rds")) %>%
    summary()
  
  # (error in here if you try to ask draws from a binomial model for rho)
  rho <- draws$quantiles[c("rho"), "50%"]
  message(paste0("Rho: ", rho))
  
  mut_data <- mut_data[!is.na(mut_data$pred),]
  
  widths <- sapply(1: floor(length(probs) / 2), 
                   function(i){probs[length(probs) - i + 1] - probs[i]})
  
  quants <- sapply(1:nrow(mut_data), function(i){
    # now assess quantiles of posterior samples, given predicted prevalences
    sims <- betabinomial_p_rho(mut_data$tested[i], mut_data$pred[i], rho) %>%
      calculate(nsim = nsim) %>%
      unlist() %>%
      quantile(probs = probs)
  }) %>%
    t()
  
  dat <- bind_cols(mut_data, quants)
  
  skip <- 7 # number of columns to skip
  ind <<- 1 
  coverages <- lapply(widths, function(width){
    if (width == 0){return(0)}
    
    lower <- skip + ind
    upper <- length(names(dat)) - ind
    # what proportion of recorded numbers of positives are in the intervals we 
    # just simulated?
    out <- sum(dat$present >= dat[,lower] & dat$present <= dat[,upper])
    
    ind <<- ind + 1
    
    out
  })
  
  data.frame(widths = widths, 
             cover = unlist(coverages) / nrow(dat))
}



posterior_predictive_check <- function(mut_data,
                                       draws_path,
                                       nsim = 500){
  # given predicted prevalences, tested, and rho, 
  # generate 100 observations from corresponding bbinomial
  # "what fraction are below the observed"?
  # histogram of the proportions - should be uniform
  
  draws <- readRDS(paste0(draws_path, "draws.rds")) %>%
    summary()
  # error in here if you try to ask draws from a binomial model for rho
  rho <- draws$quantiles[c("rho"), "50%"]
  
  message(paste0("rho ", rho))
  
  mut_data <- mut_data[!is.na(mut_data$pred),]
  
  props <- sapply(1:nrow(mut_data), function(i){
    sims <- betabinomial_p_rho(mut_data$tested[i], mut_data$pred[i], rho) %>%
      calculate(nsim = nsim)
    c(sum(sims$. < mut_data$present[i]) / nsim, 
      sum(sims$. == mut_data$present[i]) / nsim)
  })
  
  probs <- props %>%
    t() %>%
    as.data.frame() %>%
    mutate(V3 = 1 - V1 - V2) %>%
    rename(less_than = V1, equal_to = V2, greater_than = V3) %>%
    mutate(leq = less_than + equal_to, 
           geq = greater_than + equal_to,
           tested = mut_data$tested,
           present = mut_data$present,
           pred = mut_data$pred) %>%
    pivot_longer(cols = -c(tested, present, pred))
  
  return(probs)
}

# ggsave("figures/resid/residuals_k13m_bin.png", height = 3.7, width = 5, scale = 1.5)
# obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
#                "output/k13_marcse/bb_gne/preds_medians.tif",
#                xlim = c(0, 0.6), ylim = c(0, 0.6),
#                ave_tag = "_50", buffer = 100000, bb = c(27, 37, -5,  5))
# ggsave("figures/resid/residuals_k13m_bb.png", height = 3.7, width = 5, scale = 1.5)
# 
# obs_prev_panel(get_input_dir("mdr184"),
#                "output/mdr184/gneiting_sparse/preds_medians.tif",
#                xlim = c(0, 1), ylim = c(0, 1),
#                ave_tag = "_50", buffer = 100000, bb = c(27, 37, -5,  5))
# obs_prev_panel(get_input_dir("mdr184"),
#                "output/mdr184/bb_gne/preds_medians.tif",
#                xlim = c(0, 1), ylim = c(0, 1),
#                ave_tag = "_50", buffer = 100000, bb = c(27, 37, -5,  5))
# ggsave("figures/resid/residuals_mdr184_bb.png", height = 3.7, width = 5, scale = 1.5)
# 
# obs_prev_panel(get_input_dir("mdr86"),
#                "output/mdr86/gneiting_sparse/preds_medians.tif",
#                xlim = c(0, 1), ylim = c(0, 1),
#                ave_tag = "_50", buffer = 100000, bb = c(27, 37, -5,  5))
# obs_prev_panel(get_input_dir("mdr86"),
#                "output/mdr86/bb_gne/preds_medians.tif",
#                xlim = c(0, 1), ylim = c(0, 1),
#                ave_tag = "_50", buffer = 100000, bb = c(27, 37, -5,  5))
# ggsave("figures/resid/residuals_mdr86_bb.png", height = 3.7, width = 5, scale = 1.5)
# 
# mdr86

# # could try giving it longer ?
# obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
#                "output/k13_marcse/bb_gne/preds_all.tif",
#                #main = "k13 betabinom gneiting",
#                xlim = c(0, 0.6), ylim = c(0, 0.6),
#                buffer = 100000)
# ggsave("~/Desktop/presentations/MARCSE/op_bbinom.png", height=3, width=4, scale=2)
# # bit spooked by the points changing between these two ...
# # might be points falling off the mask?
# # that is so many points !
# obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
#                "output/k13_marcse/gneiting_ahmc/preds_all.tif",
#                main = "",
#                xlim = c(0, 0.6), ylim = c(0, 0.6),
#                buffer = 100000)
# ggsave("~/Desktop/presentations/MARCSE/op_binom.png", height=3, width=5, scale=1.5)



