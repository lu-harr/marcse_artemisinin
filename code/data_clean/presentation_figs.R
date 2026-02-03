# does data figs pre and post merge for astmh slides
# partner drug models are from the policy brief !
library(tidyverse)
library(terra)
library(viridisLite)
library(iddoPal)
setwd("Desktop/MARCSE/k13_seafrica")

# here's line plot over histogram of effort
markers_effort_over_time <- function(dat, 
                                     to_show = c(),
                                     bg_scale = 140,
                                     bg_col = "grey65",
                                     tag = "IDDO Surveyor"){
  markers <- dat %>%
    filter(mutant) %>% # no filter on Tested
    group_by(year, Marker) %>%
    summarise(present = sum(Present)) %>%
    rename(marker = Marker)
  
  # put test intensity in the background
  bg <- dat %>%
    group_by(Longitude, Latitude, year, PubMedID) %>%
    summarise(Tested = sum(unique(Tested))) %>%
    ungroup() %>%
    group_by(year) %>%
    summarise(Tested = sum(Tested))
  
  if (length(to_show) == 0){
    markers_keep <- markers %>%
      filter(present > 5) %>%
      arrange(present) %>%
      ungroup() %>%
      dplyr::select(marker) %>%
      unique()
  } else {
    markers_keep <- data.frame(marker = to_show)
  }
  
  
  npal = 7
  # df for line colours and styles
  sty <- markers_keep %>%
    mutate(col = factor(rep(1:npal, 5)[1:length(marker)]),
           linet = factor(rep(c(1:2), each = npal)[1:length(marker)]))
  
  # trim out markers df for line graph
  df <- markers %>%
    filter(marker %in% markers_keep$marker) %>%
    left_join(sty, join_by(marker==marker)) %>%
    mutate(group = interaction(col, linet, sep = " / "))
  
  # organise markers df for faceted plot
  # markers_keep <- markers %>%
  #   ungroup() %>%
  #   group_by(marker) %>%
  #   summarise(n = sum(present)) %>%
  #   arrange(desc(n)) %>%
  #   dplyr::slice(1:5)
  
  iddoblue <- iddoPal::iddo_palettes$iddo[1]
  
  df <- df %>%
    filter(year >= 2005) #%>%
    # left_join(markers_keep, by = join_by(marker==marker)) %>%
    # arrange(desc(n))
  
  bg <- bg %>%
    filter(year >= 2005)
  
  p1 <- 
    ggplot() +
    geom_bar(data = bg, aes(x = year, y = Tested / bg_scale), 
             stat = "identity", fill = "grey75") +
    geom_line(data = df, size = 0.7,
              aes(x = year, y = present, group = marker, color = marker, linetype = marker)) +
    geom_point(data = df, 
               aes(x = year, y = present, group = marker, color = marker)) +
    labs(
      title = paste0("Detected mutations by year - ", tag),
      color = "Marker",
      linetype = "Marker"
    ) +
    scale_color_manual(values = rep(c(iddoblue, "#c7047c", viridis(4), "#E37210"), 2)) +
    scale_linetype_manual(values = rep(1:2, each = 7)) +
    scale_y_continuous(sec.axis = sec_axis(~.*bg_scale, name="Number of tests"),
                       limits = c(0, 430)) +
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
    scale_x_continuous(breaks = seq(2006, 2024, 2)) 
  p1
}



map_k13_agg <- function(dat, tag = ""){
  dat <- dat %>%
    filter(year >= 2000) %>%
    mutate(year_bin = cut(year, breaks = seq(2000, 2025, 5))) %>%
    arrange(Present/Tested)
  
  ggplot() + 
    geom_sf(data = afr, fill = "white") + 
    geom_point(data = filter(dat, Present == 0), 
               mapping = aes(x = Longitude, y = Latitude, 
                             size = Tested, col = "grey50"),
               fill = "grey60",  pch = 21, alpha = 0.5, stroke = 0.2) +
    geom_point(data = filter(dat, Present > 0), 
               mapping = aes(x = Longitude, y = Latitude, 
                             size = Tested,
                             fill = Present / Tested),
               col = "grey50", pch=21, stroke = 0.2) +
    scale_color_manual(name = "Absence", values = c("grey30"), labels=c("")) +
    scale_fill_viridis_c(name = "Prevalence", trans = "sqrt",
                         breaks = c(0.1, 0.2, 0.4, 0.8),
                         limits = c(0, 0.81)) +
    scale_size_continuous(name = "Tested", range = c(0.5, 8), 
                          trans = "sqrt", 
                          breaks = c(10, 100, 1000, 2500), 
                          limits = c(10, 3800)) +
    #facet_wrap(~ year_bin, ncol=3) +
    labs(title = "Prevalence of Kelch 13 markers",
         subtitle = tag) +
    xlab("Longitude") +
    ylab("Latitude") +
    scale_x_continuous(breaks = seq(-20, 40, 20)) +
    scale_y_continuous(breaks = seq(-20, 40, 20)) +
    theme_grey() +
    theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(order = 2),
           size = guide_legend(order = 1))
}

moldm <- read.csv("data/clean/moldm_k13_nomarker.csv") 
map_k13_agg(moldm, tag = "IDDO Surveyor")
moldm_marcse <- read.csv("data/clean/moldm_marcse_k13_nomarker.csv")
map_k13_agg(moldm_marcse, tag = "IDDO Surveyor & MARCSE-Africa")

# I give in ... let's hard-code the markers
moldm_marcse <- read.csv("data/clean/moldm_marcse_with_markers.csv") 
moldm <- read.csv("data/clean/moldm_with_markers.csv") %>%
  filter(Tested > 5)
moldm$Marker <- moldm$strip_marker

# retrieved markers to show from moldm_marcse
markers_effort_over_time(moldm, to_show = c("C469F", "R515K", "P441L", "R539T",
                                            "R561H", "C469Y", "A675V", "M476I",
                                            "R622I", "P574L"), bg_scale = 100)
ggsave("~/Desktop/presentations/MARCSE/moldm_k13_time.png", scale = 0.8, width = 10, height = 5)

markers_effort_over_time(moldm_marcse, 
                         tag = "IDDO Surveyor & MARCSE-Africa",
                         to_show = c("C469F", "R515K", "P441L", "R539T",
                                     "R561H", "C469Y", "A675V", "M476I",
                                     "R622I", "P574L"), bg_scale = 100)
ggsave("~/Desktop/presentations/MARCSE/moldm_marcse_k13_time.png", scale = 0.8, width = 10, height = 5)


moldm <- read.csv("data/clean/moldm_k13_nomarker.csv") %>%
  filter(Present / Tested < 0.8) # suspect
# a little concerned about this picture ... 80% in Zambezi?
map_k13_agg(moldm, tag = "IDDO Surveyor")
ggsave("~/Desktop/presentations/MARCSE/moldm_k13.png", height = 6, width = 8)

moldm_marcse <- read.csv("data/clean/moldm_marcse_k13_nomarker.csv")
map_k13_agg(moldm_marcse, tag = "IDDO Surveyor & MARCSE-Africa")
ggsave("~/Desktop/presentations/MARCSE/moldm_marcse_k13.png", height = 6, width = 8)



