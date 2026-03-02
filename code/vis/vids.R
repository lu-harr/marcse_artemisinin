# make up some videos (gifs) of data/models
library(tidyverse)
library(magick)

source("code/setup.R")


make_frames <- function(marker){
  # make a map for each year we have predictions for and save em all in common dir
  
  # this was something I was sketching in to have subfigure showing data
  # if (p1 == TRUE){
  #   mut_data <- setup_mut_data(data_path, min_year = MIN_YEAR, buffer = BUFFER)
  #   
  #   p1 <- ggplot(mut_data) +
  #     geom_sf(data = afr) +
  #     geom_point(aes(x = Longitude, y = Latitude, size = Tested, col = Present / Tested)) +
  #     
  # }
  
  preds <- rast(paste0("output/", marker, "/bb_gne/preds_medians.tif"))
  pal <- ifelse(marker == "k13_marcse",
                viridis(100),
                blrd)
  
  coords <- xyFromCell(preds, cells(preds))
  vals <- terra::extract(preds, coords)
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = str_extract(lyr, "\\d{4}"))
  
  if (marker == "k13_marcse"){
    lyr_lims <- range(df$val)
    scale_trans <- "sqrt"
    scale_breaks <- c(0.1, 0.2, 0.4, 0.6)
  } else {
    lyr_lims <- c(0,1)
    scale_trans <- "identity"
    scale_breaks <- c(0, 0.25, 0.5, 0.75, 1)
  }
  
  for (yr in unique(df$year)){
    ggplot() +
      geom_sf(data = st_as_sf(afr), fill = "white") +
      geom_tile(aes(x = x, y = y, fill = val), 
                data = df %>% filter(year == yr)) +
      geom_sf(data = st_as_sf(afr), fill = NA) +
      scale_fill_viridis_c(na.value = NA, 
                           "Prevalence", 
                           trans = scale_trans, 
                           breaks = scale_breaks,
                           limits = lyr_lims) +
      scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
      scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
      labs(title = paste0(nice_name_lookup_all[marker], " prevalence - ", yr)) +
      theme_bw()
    
    ggsave(paste0("vids/", marker, "/", yr, ".png"),
           height = 6,
           width = 6)
  }
}


compile_gif <- function(marker, years = NULL, fps = 1){
  # now whack it into a gif
  pngs <- list.files(paste0("vids/", marker),
                          pattern = "\\.png$",
                          recursive = FALSE,
                          all.files = FALSE,
                          full.names = TRUE)
  
  if (!is.null(years)){
    # (implemented to chop off start of k13)
    years <- as.character(years[1]:years[2])
    pngs <- pngs[str_extract(png_files, "\\d{4}") %in% years]
  }
  
  pngs %>%
    map(image_read) %>% # reads each path file
    image_join() %>% # joins image
    image_animate(fps = fps) %>% # animates
    image_write(paste0("vids/", marker, ".gif"))
  
}


make_frames("k13_marcse")
compile_gif("k13_marcse", years = c(2010, 2028), fps = 3)


