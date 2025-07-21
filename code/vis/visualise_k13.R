# do some visualisation in here
library(viridisLite)
library(sf)
library(cowplot)

source("code/setup.R")
source("code/build_design_matrix.R")

mut_data <- setup_mut_data("data/clean/moldm_k13_nomarker.csv")
preds <- rast("output/k13/circmat_sparse/preds_all.grd")
preds <- rast("output/k13/gneiting_sparse/preds_all.grd")
test_dens <- rast("output/k13/surveillance_effort_k13.grd")

library(iddoPal)
iddoblue <- iddo_palettes_discrete$iddo[1]
# taking some bright colours that don't coincide with viridis:
case_pal <- c("#E37210", iddoblue, "#c7047c")

# mask
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

###############################################################################
# surveillance effort

# multipanel: time (hist: number tested, number of points)
surveil <- xyFromCell(test_dens, cell = cells(test_dens)) %>%
  as.data.frame() %>%
  mutate(effort = unlist(extract(test_dens, cells(test_dens))) / 4)

p1 <- ggplot(data = mut_data %>%
               group_by(year) %>%
               summarise(n = n())) +
  geom_bar(stat = "identity", aes(x = year, y = n)) +
  ylab("Number of locations") +
  xlab("Year") +
  ggtitle("(a)")

p2 <- ggplot(data = mut_data %>%
               group_by(year) %>%
               summarise(tested = sum(tested))) +
  geom_bar(stat = "identity", aes(x = year, y = tested)) +
  ylab("Number of tests") +
  xlab("Year") +
  ggtitle("(b)")

# could just go back to totals? As emphasis is on effort?
# p2 <- ggplot(data = mut_data %>%
#                group_by(year) %>%
#                summarise(present = sum(present),
#                          absent = sum(tested) - present) %>%
#                pivot_longer(!year, names_to = "Tests", values_to = "Count")) +
#   geom_bar(stat = "identity", aes(x = year, y = Count, fill = Tests)) +
#   xlab("Year")

p3 <- ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") + # not showing anything in the background here ...
  geom_tile(data = surveil, aes(x, y, fill = effort)) +
  geom_sf(data = st_as_sf(afr), colour = "white", fill = NA) +
  scale_fill_viridis_c(na.value = NA, bquote(atop("Tests per","~100"~km^2)), trans="sqrt") +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("(c)")

library(deeptime)

gg1 <- ggarrange2(p1, p2, layout = rbind(c(1), c(2)), draw = FALSE)
ggarrange2(gg1, p3, widths = c(1,2))
ggsave("figures/surveillance_effort_k13.png", ggarrange2(gg1, p3, widths = c(1,2)),
       height = 5, width = 7.5, scale = 1.5)
# don't love that I'm judging vertical stretching by eye
# It only took me 45 mins to work this out I guess

###############################################################################
# visualise lengthscale priors:

# axis(xxx, labels = degrees_to_radians(xxx), "Distance (degrees)")
# axis(xxx, labels = distHaversine(xxx), "Distance (km)")
library(geosphere)
library(greta.gp)

circmat_len <- lognormal(meanlog = 0, sdlog = 1)
xxx = seq(0, 1, length.out = 100)
sdlogs = c(2.5, 2, 1.5, 1)
sdlogs = rep(1, 3)
mulogs = c(-1,-0.5, 0, 0.5)

dist_rads <- seq(0, 0.5, length.out=6)
dist_degs <- radians_to_degrees(dist_rads)
dist_ms <- distHaversine(c(0,0), matrix(c(rep(0, 6), dist_degs), ncol = 2))

cand_rads <- c(0, 0.0785, 0.158, 0.235, 0.314, 0.392)#, 0.472)
cand_degs <- radians_to_degrees(cand_rads)
dist_ms <- distHaversine(c(0,0), matrix(c(rep(0, 6), cand_degs), ncol = 2))
dist_ms

draws <- read_rds("output/circmat_k13/draws.rds")
summ_s <- summary(draws$`11`[,1])
summ_t <- summary(draws$`11`[,3])

par(mfrow = c(1,2), oma = c(8,0,0,0))

prior_meanlog = -2
prior_sdlog = 1

scaled_years <- scale_years(range(pfpr_years)) %>%
  unlist()

# can't really think any more tbh
plot(xxx, dlnorm(xxx, meanlog = prior_meanlog, sdlog = prior_sdlog), 
     xlab = "Lengthscale (scaled_years)", ylab = "p(lengthscale)", 
     main = "Lognormal prior on temporal lengthscale", type="l", xlim = c(0, 5))
axis(1, 0:5 * 0.147442, 0:5, line = 5)

plot(xxx, dlnorm(xxx, meanlog = prior_meanlog, sdlog = prior_sdlog), 
     xlab = "Lengthscale (radians)", ylab = "p(lengthscale)", 
     main = "Lognormal prior on spatial lengthscale", type="l", xlim = c(0, 0.42))
#axis(1, dist_rads, labels = round(dist_degs, digits=2), line = 5)
axis(1, degrees_to_radians(seq(0, 25, length.out = 6)), 
     labels = seq(0, 25, length.out = 6), line = 5)
mtext(side = 1, "Lengthscale (degrees)", line = 8)
axis(1, cand_rads, labels = seq(0, 2500, length.out = 6), #round(dist_ms/1000), 
     line = 10)
mtext(side = 1, "Lengthscale (km)", line = 13)
abline(v = summ_s$quantiles[3])
abline(v = summ_s$quantiles[c(1,5)], lty=2)
legend("topright", lty=1:2, c("Median", "2.5% - 97.5%"), title = "Draws")
# add priors to hists ...?

###############################################################################
# require common colour palette between years !

zambezi <- list(xmin = 19, xmax = 26, ymin = -21, ymax = -14)
victoria <- list(xmin = 27, xmax = 35, ymin = -6, ymax = 2)
eswatini <- list(xmin = 34, xmax = 41, ymin = 12, ymax = 19)

zoom_df <- rbind(zambezi, victoria, eswatini) %>%
  as.data.frame() %>%
  unnest() %>%
  suppressWarnings()

gg_ras_prep <- function(ras, extent = NULL, shp = NULL){
  # this was annoying me
  if (!is.null(extent)){
    extent <- ext(unlist(extent))
    ras <- mask(ras, extent) %>%
      crop(extent)
    
    if (!is.null(shp)){
      # assuming we're working with sf
      shp <- st_crop(shp, extent)
    }
  }
  
  coords <- xyFromCell(ras, cells(ras))
  vals <- terra::extract(ras, coords)
  
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4),
           tag = substr(lyr, 11, 14))
  
  list(df = df,
       shp = shp)
}

df <- gg_ras_prep(preds)$df
zambezi_bits <- gg_ras_prep(preds, zambezi, afr)
victoria_bits <- gg_ras_prep(preds, victoria, afr)
eswatini_bits <- gg_ras_prep(preds, eswatini, afr)

years_to_plot <- c("2014", "2018", "2022")

zoom_pan <- function(bits, 
                     yearr, 
                     tagg = "medi",
                     pal = viridis(100), 
                     scale_lims = c(0, 0.64), #c(0, 0.42),
                     panel_col = "black"){
  message(yearr)
  ggplot() +
    geom_tile(data = bits$df %>%
                filter(year == yearr & tag == tagg), 
              mapping = aes(x = x, y = y, fill = val)) +
    geom_sf(data = bits$shp, fill = NA, col = "grey") +
    scale_fill_gradientn(colours = pal, limits = scale_lims, trans = "sqrt") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          legend.justification = "top",
          axis.title = element_blank(),
          legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(rep(0.2,4), "cm"),
          panel.border = element_rect(colour = panel_col,
                                      linewidth = 1.5))
}


medians <- ggplot() +
  geom_sf(data = afr, fill = "white") +
  geom_tile(data = df %>%
              filter(year %in% years_to_plot & tag == "medi"),
            mapping = aes(x = x, y = y, fill = val)) +
  facet_wrap(~year, ncol = 1, strip.position = "left") +
  # geom_tile(data = df %>%
  #             filter(year =="2022" & tag == "medi"),
  #           mapping = aes(x = x, y = y, fill = val)) +
  scale_fill_viridis_c(na.value = NA, "Prevalence", trans = "sqrt") +
  # xlab("Longitude") +
  # ylab("Latitude") +
  # labs(title = "Median") +
  geom_rect(zoom_df %>% mutate(year = years_to_plot[1]), 
            mapping = aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax), colour = case_pal, linewidth = 1, fill = NA) +
  geom_text(data = data.frame(lab = c("a", "b", "c"), 
                              year = c("2014", "2018", "2022")), 
            mapping = aes(x = 0, y = 0, label = lab)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(), #element_text(hjust = 0.5),
        legend.justification = "top") 
medians

# yr1 <- plot_grid(zoom_pan(zambezi_bits, "2014", panel_col = case_pal[1]),
#                  zoom_pan(victoria_bits, "2014", panel_col = case_pal[2]),
#                  NULL,
#                  zoom_pan(eswatini_bits, "2014", panel_col = case_pal[3])) +
#   theme(plot.margin = unit(rep(0.2,4), "cm"))
# 
# yr2 <- plot_grid(zoom_pan(zambezi_bits, "2018", panel_col = case_pal[1]),
#                  zoom_pan(victoria_bits, "2018", panel_col = case_pal[2]),
#                  NULL,
#                  zoom_pan(eswatini_bits, "2018", panel_col = case_pal[3])) +
#   theme(plot.margin = unit(rep(0.2,4), "cm"))
# 
# yr3 <- plot_grid(zoom_pan(zambezi_bits, "2022", panel_col = case_pal[1]),
#                  zoom_pan(victoria_bits, "2022", panel_col = case_pal[2]),
#                  NULL,
#                  zoom_pan(eswatini_bits, "2022", panel_col = case_pal[3])) +
#   theme(plot.margin = unit(rep(0.2,4), "cm"))
# zooms <- plot_grid(yr1, yr2, yr3, ncol = 1)


zooms <- plot_grid(zoom_pan(eswatini_bits, "2014", panel_col = case_pal[3]),
                   zoom_pan(victoria_bits, "2014", panel_col = case_pal[2]),
                   zoom_pan(zambezi_bits, "2014", panel_col = case_pal[1]),
                   zoom_pan(eswatini_bits, "2018", panel_col = case_pal[3]),
                   zoom_pan(victoria_bits, "2018", panel_col = case_pal[2]),
                   zoom_pan(zambezi_bits, "2018", panel_col = case_pal[1]),
                   zoom_pan(eswatini_bits, "2022", panel_col = case_pal[3]),
                   zoom_pan(victoria_bits, "2022", panel_col = case_pal[2]),
                   zoom_pan(zambezi_bits, "2022", panel_col = case_pal[1]),
                   ncol = 1) +
  theme(plot.margin = unit(rep(0.1,4), "cm"))



sds <- ggplot() +
  geom_sf(data = afr, fill = "white") +
  geom_tile(data = df %>%
              filter(year %in% years_to_plot & tag == "sd"), 
            mapping = aes(x = x, y = y, fill = val)) +
  facet_wrap(~year, ncol = 1) +
  scale_fill_distiller(palette = "Oranges", 
                       na.value = NA, 
                       "Uncertainty", 
                       direction = 1,
                       trans = "sqrt") +
  # xlab("Longitude") +
  # ylab("") +
  # labs(title = "Standard deviation") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_blank(), # element_text(hjust = 0.5),
        legend.justification = "top") 

legs <- plot_grid(get_legend(medians),
                  get_legend(sds),
                  # this is a bit hacky
                  NULL,
                  NULL,
                  ncol = 1)

df <- data.frame(xmin = rep(0.001, 3),
                 xmax = rep(0.028, 3),
                 ymin = c(0.0155, 0.337, 0.663),
                 ymax = c(0.315, 0.639, 0.965),
                 lab = c("D1246Y", "Y184F", "N86Y"))
p <- p + 
  # tried adding outer margin but that didn't do anything
  geom_rect(data = df, aes(xmin=xmin, xmax=xmax, ymin=ymin, 
                           ymax=ymax), 
            colour="grey10", fill="grey85", linewidth=0.3) +
  geom_text(data = df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                           label = lab), angle = 90)
  
df <- data.frame(xmin = c(0.04, 0.525),
                 xmax = c(0.5, 0.873),
                 ymin = rep(0.992, 2),
                 ymax = rep(1.02, 1),
                 lab = c("Median", "Standard deviation"))

plot_grid(medians + theme(legend.position = "none"), 
          zooms, 
          sds + theme(legend.position = "none"), 
          legs,
          ncol = 4, rel_widths = c(1,0.3,0.933,0.3)) +
  theme(plot.margin = unit(c(0.7,0,0,0), "cm")) +
  geom_rect(data = df, aes(xmin=xmin, xmax=xmax, ymin=ymin, 
                           ymax=ymax), 
            colour="grey10", fill="grey85", linewidth=0.3) +
  geom_text(data = df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                           label = lab))

ggsave("figures/k13_out_gneiting.png", height = 5.2, width = 4.5, scale = 2)



################################################################################
# WITHOUT ZOOM PANS

medians <- ggplot() +
  geom_sf(data = afr, fill = "white") +
  geom_tile(data = df %>%
              filter(year %in% years_to_plot & tag == "medi"),
            mapping = aes(x = x, y = y, fill = val)) +
  facet_wrap(~year, ncol = 1, strip.position = "left") +
  scale_fill_viridis_c(na.value = NA, "Prevalence") +theme_bw() +
  labs(title = "Median") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #plot.title = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.justification = "top")

medians

sds <- ggplot() +
  geom_sf(data = afr, fill = "white") +
  geom_tile(data = df %>%
              filter(year %in% years_to_plot & tag == "sd"), 
            mapping = aes(x = x, y = y, fill = val)) +
  facet_wrap(~year, ncol = 1) +
  scale_fill_distiller(palette = "Oranges", 
                       na.value = NA, 
                       "Uncertainty", 
                       direction = 1,
                       trans = "sqrt") +
  labs(title = "Standard deviation") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        #plot.title = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.justification = "top") 
sds


