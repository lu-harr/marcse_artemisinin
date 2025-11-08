library(viridisLite)
library(RColorBrewer)
library(iddoPal)
source("code/setup.R")
source("code/vis/vis_funcs.R")

oranges <- brewer.pal(9, "Oranges")
blrd <- iddoPal::iddo_palettes$BlGyRd
iddoblue <- iddo_palettes_discrete$iddo[1]
# taking some bright colours that don't coincide with viridis:
case_pal <- c("#E37210", iddoblue, "#c7047c")

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
                             exte,
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
      filter(year > 2005 & present > 0) %>%
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
    filter(year > 2005 & present > 0) %>%
    group_by(marker) %>%
    summarise(n = sum(present)) %>%
    arrange(desc(n))
  
  nnonzero <- nrow(markers_keep)
  markers_keep <- filter(markers_keep, n > 0) %>%
    dplyr::slice(1:npal)
  
  if (nnonzero > npal){
    markers_keep <- dplyr::slice(markers_keep, 1:(npal - 1)) %>%
      bind_rows(data.frame(marker = "Others", n=1))
  }
    
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
    facet_wrap(vars(Marker), ncol = 4) +
    # labs(title = "(b) Prevalence of k13 markers in Africa") +
    xlab("Longitude") +
    ylab("Latitude") +
  # scale_x_continuous(breaks = seq(-20, 40, 20)) +
  # scale_y_continuous(breaks = seq(-20, 40, 20))
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing = unit(0, "lines"))
  
  
  bg_col <- "grey65"
  present_lims <- c(0, max(df$present) %/% 10 * 10 + 10)
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


# perhaps this should be wrapped into its own script ...
# make calls to ribbon plot, row plot, etc., but provide masked raster

country_profile <- function(iso = c("KEN"), # of values in afr$iso_a3; could be vector for regional map
                            mod = c("gneiting_sparse", "bb_gne"), 
                            marker = c("k13_marcse", "partner"),
                            buff = NULL, # option to look at neighbouring countries; in m
                            epsg = 32736,
                            years = as.character(seq(2009, 2024, length.out = 4)),
                            summarise_only = FALSE){
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
    paste0("/preds_medians.tif") %>%
    rast() %>%
    #aggregate(fact = 10) %>%
    mask(vect(buff)) %>%
    trim()
  
  if (summarise_only){
    preds <- preds[[str_extract(names(preds), "\\d{4}") %in% years]]
    
    # plot(preds)
    # mtext(marker)
    
    medians <- terra::global(preds, fun = "mean", na.rm = TRUE)
    
    return(data.frame(year = years,
                      median = medians))
  }
  
  # top panel: data
  # if (marker == "k13_marcse"){
    dat_plot <- plot_k13_markers(buff, exte, markers_keep = NULL)
    pred_plot <- map_pred_row(get_output_dir(marker, mod) %>%
                                paste0("/preds_medians.tif"),
                              field = "50",
                              pal = viridis(10),
                              buff = buff,
                              exte = exte,
                              years = years, 
                              top_pan = TRUE, 
                              ylab = "")
  # } else {
  #   # partner drugs data
  #   dat_plot <- plot_partner_markers(buff, marker)
  #   # need to wrap into markers?
  #   pred_plot <- map_pred_row()
  # }
  
  list(dat_plot, pred_plot)
}

# out = country_profile(iso = c("RWA", "UGA"),
#                       marker = "k13_marcse", 
#                       mod = "bb_gne", 
#                       buff = 50000)
# 
# p <- plot_grid(out[[1]][[2]], out[[1]][[1]], out[[2]], 
#                ncol = 1, rel_heights = c(0.4, 1, 0.55))
# ggsave("figures/country_profiles/rwa_uga.png", p, height = 10, width = 9)



# out = country_profile(iso = c("RWA", "UGA", "KEN"),
#                       marker = "k13_marcse", 
#                       mod = "bb_gne", 
#                       buff = 50000)
# 
# p <- plot_grid(out[[1]][[2]], out[[1]][[1]], out[[2]], 
#                ncol = 1, rel_heights = c(0.7, 1, 0.55))
# ggsave("figures/country_profiles/rwa_uga_ken.png", p, height = 9, width = 9)



# out = country_profile(iso = c("RWA", "UGA", "BDI"),
#                       marker = "k13_marcse", 
#                       mod = "bb_gne", 
#                       buff = 50000)
# 
# p <- plot_grid(out[[1]][[2]], out[[1]][[1]], out[[2]], 
#                ncol = 1, rel_heights = c(0.48, 1, 0.5))
# ggsave("figures/country_profiles/rwa_uga_bdi.png", p, height = 12, width = 9)


# out = country_profile(iso = c("ETH", "ERI", "SSD", "SDN"),
#                       marker = "k13_marcse", 
#                       mod = "bb_gne", 
#                       buff = 50000)
# 
# p <- plot_grid(out[[1]][[2]], out[[1]][[1]], out[[2]], 
#                ncol = 1, rel_heights = c(0.48, 1, 0.5))
# ggsave("figures/country_profiles/northern_cluster.png", p, height = 7.5, width = 9)





