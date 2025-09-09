# folding functions for plotting in here
library(looseVis)
library(cowplot)
library(iddoPal)
library(sf)

# wrapping up plot of all-Africa posterior medians over times
pred_time_plot <- function(in_path, 
                           title = "",
                           points_path = "",
                           pal = iddoPal::iddo_palettes$soft_blues,
                           zooms = NULL,
                           zoom_pal = NULL, 
                           alpha = 0.5){
  
  preds <- rast(in_path) %>%
    aggregate(fact = 10)
  message("aggregating")
  
  coords <- xyFromCell(preds, cells(preds))
  vals <- terra::extract(preds, coords)
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4),
           tag = substr(lyr, 6, 14)) %>% # pick out year and thingo
    filter(tag == "50") %>%
    group_by(year) %>%
    summarise(q = list(quantile(val, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)))) %>%
    unnest_wider(q) %>%
    ungroup() %>%
    mutate(year = as.numeric(year))
  
  p <- ggplot(df) +
    geom_line(aes(x = year, y = `0%`, linetype = "0% - 100%")) +
    geom_line(aes(x = year, y = `100%`, linetype = "0% - 100%")) +
    geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = "2.5% - 97.5%"), alpha = alpha) + #fill=pal[6]) +
    geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`, fill = "25% - 75%"), alpha = alpha) + #fill=pal[4]) +
    geom_ribbon(aes(x = year, ymin = `50%`, ymax = `50%`, fill = "50%")) + #fill=pal[1]) +
    geom_line(aes(x = year, y = `50%`), col = pal[1], linewidth = 1) +
    scale_linetype_manual("", values = c("0% - 100%" = 2)) +
    scale_fill_manual("", values = c("2.5% - 97.5%" = pal[6], "25% - 75%" = pal[4], "50%" = pal[1])) +
    ylab("Prevalence") +
    xlab("Year") +
    labs(title = title) +
    ylim(0, 1) +
    scale_x_continuous(breaks = seq(2000, 2022, 2), expand = c(0,0)) +
    theme_bw() +
    theme(legend.spacing.y = unit(-10, "cm"),
          legend.background = element_rect(fill = NA))
  
  if (!is.null(zooms)){
    for (i in 1:nrow(zooms)){
      exte <- ext(unlist(zooms[i,]))
      tmp <- crop(preds, exte)
      coords <- xyFromCell(tmp, cells(tmp))
      vals <- terra::extract(tmp, coords)
      
      df <- cbind(coords, vals) %>%
        pivot_longer(starts_with("2"),
                     names_to = "lyr",
                     values_to = "val") %>%
        mutate(year = substr(lyr, 1, 4),
               tag = substr(lyr, 11, 14)) %>% # pick out year and thingo
        filter(tag == "medi") %>%
        group_by(year) %>%
        summarise(med = median(val)) %>%
        ungroup() %>%
        mutate(year = as.numeric(year))
      
      p <- p + geom_line(df, mapping = aes(x = year, y = med), colour = case_pal[i])
    }
  }
  
  
  
  if (points_path != ""){
    message("Watch out! I set limits manually!")
    mut_data <- setup_mut_data(points_path, min_year = MIN_YEAR)
    p <- p + geom_point(aes(x=jitter(year), y=present/tested,  
                            size=tested), 
                        colour="grey", pch = 21,
                        mut_data) +
      scale_size_continuous(name = "Tested", trans = "sqrt", 
                            range = c(0.2, 4), limits = c(5, 3500)) # +
    # geom_boxplot(aes(x = year, y = present/tested, group = as.factor(year)),
    #              mut_data, fill = NA, outliers = FALSE)
  }
  
  p
}


# to wrap up prediction figures for partner drug models :
map_pred_row <- function(in_path,
                         years,
                         pal,
                         field = c("50", "sd"),
                         buff = NULL,
                         xlab = "Longitude",
                         ylab = "Latitude",
                         legend_lim = waiver(),
                         top_pan = FALSE){
  preds <- rast(in_path)
  preds <- preds[[grep(paste0(field, "$"), names(preds))]]
  
  if (!is.null(buff)){
    preds <- mask(preds, buff) %>%
      trim()
  }
  
  coords <- xyFromCell(preds, cells(preds))
  vals <- terra::extract(preds, coords)
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4)) %>%
    filter(year %in% years) # pick out year
  
  p <- ggplot() +
    geom_sf(data = afr, fill = "white") +
    geom_tile(data = df, 
              mapping = aes(x = x, y = y, fill = val)) +
    facet_wrap(~year, nrow = 1) +
    #xlab(xlab) +
    ylab(ylab) +
    #scale_x_continuous(breaks = seq(-20, 40, 20)) +
    #scale_y_continuous(breaks = seq(-20, 40, 20)) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_blank(),
          legend.justification = "left",
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 0, hjust = 1),
          #axis.title.y = element_blank(),
          title = element_blank(),
          panel.spacing = unit(0, "lines"))
  
  if (!top_pan){
    p <- p + theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
  }
  
  if (field == "50"){
    p <- p + scale_fill_gradientn(name = "Prevalence",
                                  colors = pal, 
                                  breaks = c(0, 0.5, 1), 
                                  labels = c("0  (all wildtype)", "0.5", "1  (all mutant)"),
                                  limits = legend_lim)
  } else {
    p <- p + scale_fill_gradientn(colors = pal, 
                                  name = "Estimate SD",
                                  limits = legend_lim)
  }
  
  p
}


# observed vs predicted values
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
  mut_data <- setup_mut_data(data_path, min_year = MIN_YEAR)
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
  
  message(paste0("Sites with preds: ", nrow(mut_data)))
  message(paste0("Goodness of fit: ", 
                 sum((mut_data$pred*mut_data$tested - mut_data$present)**2 / 
                          (mut_data$pred*mut_data$tested))))
  
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
    xlab("Year") +
    ylab("Observed - Predicted prevalence") +
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
    plot_grid(p3, rel_widths = c(0.4, 0.75))
}



# # e.g.:
# # might want to re-land some points inside of model fitting
# obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
#                "output/k13_marcse/gneiting_sparse/preds_all.grd",
#                #main = "k13 gneiting binom", 
#                xlim = c(0, 0.6), ylim = c(0, 0.6),
#                ave_tag = "_post_median", buffer = 100000)
# ggsave("~/Desktop/presentations/marcse/residuals_k13m_t.png", height = 3.7, width = 5, scale = 1.5)
# 
# # could try giving it longer ?
obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
               "output/k13_marcse/bb_gne/preds_all.tif",
               #main = "k13 betabinom gneiting",
               xlim = c(0, 0.6), ylim = c(0, 0.6),
               buffer = 100000)
ggsave("~/Desktop/presentations/MARCSE/op_bbinom.png", height=3, width=4, scale=2)
# # bit spooked by the points changing between these two ...
# # might be points falling off the mask?
# # that is so many points !
obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
               "output/k13_marcse/gneiting_ahmc/preds_all.tif",
               main = "",
               xlim = c(0, 0.6), ylim = c(0, 0.6),
               buffer = 100000)
ggsave("~/Desktop/presentations/MARCSE/op_binom.png", height=3, width=5, scale=1.5)



## country-level plots
# give me directory where all of the model objects/predictions are ...
get_output_dir <- function(marker, mod){
  paste0("output/", marker, "/", mod)
}


# give me location of data to bring in ...
get_input_dir <- function(snp){
  in_dat <- ifelse(snp == "k13",
                   "data/clean/moldm_k13_nomarker.csv",
                   ifelse(snp == "crt76",
                          "data/clean/moldm_crt76.csv",
                          ifelse(snp == "k13_marcse",
                                 "data/clean/moldm_marcse_k13_nomarker.csv",
                                 paste0(paste0("data/clean/pfmdr_single_", snp, ".csv")))))
  in_dat
  
}


# for k13:
plot_k13_markers <- function(buff, 
                             markers_keep = NULL, 
                             npal = 7){
  # relies on marcse_merge - need to clean up data workflow at some point ..
  moldm <- read.csv("data/clean/moldm_marcse_with_markers.csv") %>%
    drop_na(Longitude, Latitude) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = st_crs(buff)) %>%
    st_intersection(st_geometry(buff)) %>%
    suppressWarnings()
  
  # probably a simpler way to do this ....
  moldm <- cbind(st_coordinates(moldm),
                 st_drop_geometry(moldm))
  moldm <- moldm[, !duplicated(names(moldm))]
    
  markers <- moldm %>%
    filter(mutant) %>% # no filter on Tested
    group_by(year, Marker) %>%
    summarise(present = sum(Present)) %>%
    rename(marker = Marker)
  
  wildtypes_to_add <- anti_join(
    # locations in `wildtypes` that do not occur in `mutants`,
    # paying attention to `Tested` but NOT to `pubs`
    moldm %>%
      filter(Marker == "wildtype") %>%
      dplyr::select(X, Y, year, Tested, Site.Name, Country) %>%
      unique(),
    moldm %>%
      filter(mutant) %>%
      dplyr::select(X, Y, year, Tested, Site.Name, Country)) %>%
    mutate(Present = 0) %>%
    suppressMessages()
  
  # put surveillance intensity in the background of time figure
  bg <- moldm %>%
    group_by(X, Y, year, PubMedID) %>%
    summarise(Tested = sum(unique(Tested))) %>%
    ungroup() %>%
    group_by(year) %>%
    summarise(Tested = sum(Tested)) %>%
    filter(year > 2005)
  
  # I cannot state again how intensely easy this is in base graphics
  if (is.null(markers_keep)){
    markers_keep <- markers %>%
      filter(present > 2) %>%
      arrange(present) %>%
      ungroup() %>%
      dplyr::select(marker) %>%
      unique()
  } else {
    markers_keep <- data.frame(marker = markers_keep)
  }
  
  sty <- markers_keep %>%
    mutate(col = factor(rep(1:npal, 5)[1:length(marker)]),
           linet = factor(rep(c(1:2), each = npal)[1:length(marker)]))
  
  df <- markers %>%
    filter(marker %in% markers_keep$marker) %>%
    left_join(sty, join_by(marker==marker)) %>%
    mutate(group = interaction(col, linet, sep = " / ")) %>%
    filter(year > 2005)
  
  # include slice end as func par?
  markers_keep <- markers %>% 
    ungroup() %>%
    group_by(marker) %>%
    summarise(n = sum(present)) %>%
    arrange(desc(n)) %>%
    slice(1:5) %>%
    bind_rows(data.frame(marker = "Others", n=1))
  
  wts <- wildtypes_to_add %>%
    #rename(Longitude = X, Latitude = Y) %>%
    dplyr::select(X, Y, year, Tested) %>%
    filter(year > 2008) %>%
    merge(markers_keep, by = NULL) %>%
    #rename(y = marker) %>% # what was I doing with this again?
    rename(Marker = marker) %>%
    mutate(Present = 0)
  
  markers_disagg <- moldm %>%
    filter(mutant) %>%
    mutate(Marker = ifelse(Marker %in% markers_keep$marker,
                           Marker, "Others")) %>%
    filter(!(Marker == "Others" & Present == 0) & # don't want zeroes for all other markers, just WTs as BG
             Tested > 5) %>% # don't want prevalence == 1, tested < 5 points
    bind_rows(wts) %>% #  %>% mutate(Marker = "wt")
    # add back in later?
    # mutate(year_bin = cut(year, breaks = c(min(year)-1, 2012, 2015, 2018, 2021, max(year)))) %>%
    mutate(Marker = factor(Marker, levels = markers_keep$marker))
  
  p1 <- ggplot() +
    geom_sf(data = exte, fill = "white") + 
    geom_point(data = filter(markers_disagg, Present == 0),
               mapping = aes(x = X, y = Y,
                             size = Tested, col = "grey50"),
               fill = "grey70",  pch = 21, alpha = 0.5, stroke = 0.2) +
    geom_point(data = filter(markers_disagg, Present > 0), 
               mapping = aes(x = X, y = Y, 
                             size = Tested,
                             fill = Present / Tested),
               col = "grey50", pch=21, stroke = 0.2) +
    scale_color_manual(name = "", values = c("grey30"), labels=c("Absence")) +
    scale_fill_viridis_c(name = "Prevalence", trans = "sqrt") +
    scale_size_continuous(name = "Tested", range = c(0.2, 4), trans = "sqrt") +
    facet_wrap(vars(Marker)) +
    # labs(title = "(b) Prevalence of k13 markers in Africa") +
    xlab("Longitude") +
    ylab("Latitude") #+
    # scale_x_continuous(breaks = seq(-20, 40, 20)) +
    # scale_y_continuous(breaks = seq(-20, 40, 20))
  
  
  bg_col <- "grey65"
  present_lims <- c(0, max(df$present) %/% 50 * 50 + 50)
  bg_scale <- max(bg$Tested) / max(df$present)
  
  p2 <- ggplot() +
    geom_bar(data = bg, aes(x = year, y = Tested / bg_scale), 
             stat = "identity", fill = "grey75") +
    geom_line(data = df, size = 0.7,
              aes(x = year, y = present, group = marker, color = marker, linetype = marker)) +
    geom_point(data = df, 
               aes(x = year, y = present, group = marker, color = marker)) +
    labs(
      # title = "(a) Detected mutations by year",
      color = "Marker",
      linetype = "Marker"
    ) +
    scale_color_manual(values = rep(c(viridis(4), "#E37210", iddoblue, "#c7047c"), 2)) +
    scale_linetype_manual(values = rep(1:2, each = 7)) +
    scale_y_continuous(sec.axis = sec_axis(~.*bg_scale, name="Number of tests"),
                       limits = present_lims) +
    theme_minimal() +
    xlab("Year") +
    ylab("Mutations detected") +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.key.spacing.y = unit(-0.3, "lines"),
          legend.box.margin = margin(50, 6, 6, 6),
          axis.text.y.right = element_text(color = bg_col),
          axis.title.y.right = element_text(color = bg_col),
          axis.ticks.y.right = element_line(color = bg_col)) +
    #legend.justification.right = "bottom") +
    scale_x_continuous(breaks = 2005:2024) 
  
  # could grid things at this point?
  return(list(p1, p2))
}


# unfinished: ...
plot_partner_markers <- function(shp, buff, 
                                 markers = c("mdr86", "mdr184", "mdr1246", "crt76")){
  for (marker in markers){
    mut_data <- setup_mut_data(get_input_dir(marker), min_year = MIN_YEAR)
  }
  
}

# perhaps this should be wrapped into its own script ...
# make calls to ribbon plot, row plot, etc., but provide masked raster
country_profile <- function(iso = c("KEN"), # of values in afr$iso_a3
                            # could be vector? for regional map?
                            mod = c("gneiting_sparse", "bb_gne"), 
                            marker = c("k13_marcse", "partner"),
                            buff = NULL, # option to look at neighbouring countries; in km
                            epsg = 32736){
  # a function to give me country-level plots for a given model/marker?
  # would like:
    # - map of data, hist of survey effort, number of contributing studies, number of testees
    # - map of model outputs for a couple of time-points, ribbon of medians
    # - residuals?
  
  exte <- afr %>% 
    filter(iso_a3 %in% iso) %>%
    st_transform(epsg)
  if (!is.null(buff)){
    buff <- st_buffer(exte, buff)
      # don't do this: introduces "duplicate edges"
      #st_union() # removes internal borders if we're asking for multiple countries
    # could also simplify here ...
  } else {
    buff = exte
  }
  
  buff <- st_transform(buff, 4326) # back for latlons
  exte <- st_transform(exte, 4326)
  
  preds <- get_output_dir(marker, mod) %>%
    paste0("/preds_all.tif") %>%
    rast() %>%
    aggregate(fact = 10) %>%
    mask(vect(buff)) %>%
    trim()
  
  # top panel: data
  if (marker == "k13_marcse"){
    dat_plot <- plot_k13_markers(buff)
    pred_plot <- map_pred_row(get_output_dir(marker, mod) %>%
                                paste0("/preds_all.tif"),
                              field = "50",
                              pal = viridis(10),
                              years)
  } else {
    # partner drugs data
    dat_plot <- plot_partner_markers(buff, marker)
    # need to wrap into markers?
    pred_plot <- map_pred_row()
  }
  
  
}

# country_profile(iso = c("KEN", "RWA", "UGA"),
#                 marker = "k13_marcse", mod = "bb_gne")
# 
# 
# exte = afr %>% 
#   filter(iso_a3 %in% c("KEN", "RWA", "UGA")) %>%
#   st_transform(32736)
# plot(st_geometry(exte))
# buff = st_buffer(exte, 100000)
# plot(st_geometry(buff), add=TRUE)





# this absolutely came from chatty g
# | Region                   | UTM Zone | EPSG Code   |
#   | ------------------------ | -------- | ----------- |
#   | Kenya/Tanzania           | 36S      | 32736       |
#   | Uganda/South Sudan       | 36N      | 32636       |
#   | Ethiopia                 | 37N      | 32637       |
#   | Zambia/Malawi/Mozambique | 36S–37S  | 32736–32737 |
#   | South Africa             | 34S–36S  | 32734–32736 |




