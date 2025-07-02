# do some visualisation in here
library(viridisLite)
library(sf)

source("code/setup.R")
source("code/build_design_matrix.R")

mut_data <- setup_mut_data("data/moldm_k13_nomarker.csv")
preds <- rast("output/circmat_k13/preds_all.grd")
test_dens <- rast("output/surveillance_effort_k13.grd")

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

axis(xxx, labels = degrees_to_radians(xxx), "Distance (degrees)")
axis(xxx, labels = distHaversine(xxx), "Distance (km)")

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
# nice clean traceplot
library(GGally)
library(brms)
library(bayesplot)

draws <- read_rds("output/circmat_k13/draws.rds")

post <- as_draws_df(draws) %>%
  rename("Lengthscale (spatial)" = "circmat_len",
         "Variance (spatial)" = "circmat_var", 
         "Lengthscale (temporal)" = "expo_len",    
         "Variance (temporal)" = "expo_var", # (years are scaled - unscale relevant params?)
         "Nugget variance" = "nugget_var",
         "Beta (intercept)" = "beta[1,1]",
         "Beta (scaled year)" = "beta[2,1]",
         "Beta (PfPR)" = "beta[3,1]")

color_scheme_set("purple") # the bayesplot scheme is much better suited to this application ..
bayesplot::mcmc_trace(post)
ggsave("figures/trace_k13.png", height=10, width=15)



modified_density = function(data, mapping, ...) {
  ggally_densityDiag(data, mapping, ...) + 
    scale_fill_manual(values = c("#e5cce5","#bf7fbf","#a64ca6","#800080","#660066","#400040"))
}

# ggplot is the silliest darn software going
# why would anyone want or need to change a colour palette?
# beats me
modified_points = function(data, mapping, ...) {
  ggally_points(data, mapping, ...) + 
    scale_color_manual(values = c("#e5cce5","#bf7fbf","#a64ca6","#800080","#660066","#400040"))
}

modified_cor = function(data, mapping, ...){
  ggally_cor(data, mapping, ...) +
    scale_color_manual(values = c("#e5cce5","#bf7fbf","#a64ca6","#800080","#660066","#400040"))
}

# love the purps but they're a bit too light to use for cor
post %>%
  mutate(chain = as.factor(.chain)) %>%
  ggpairs(columns = 1:6,
          mapping = aes(colour = chain),
          lower = list(continuous = wrap(modified_points, alpha = 0.4)), # can't seem to turn alpha down here?
          diag = list(continuous = wrap(modified_density, alpha = 0.5)),
          upper = list(continuous = modified_cor)) # probably don't need corrs for chains?
ggsave("figures/chain_corr_k13.png", height=12, width=12)

###############################################################################
preds <- rast("output/circmat_k13/preds_all.grd")
# require common colour palette between years !


coords <- xyFromCell(preds, cells(preds))
vals <- terra::extract(preds, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         tag = substr(lyr, 11, 14)) # pick out year and thingo

years_to_plot <- c("2014", "2018", "2022")

p1 <- ggplot() +
  geom_sf(data = afr, fill = "white") +
  geom_tile(data = df %>%
              filter(year %in% years_to_plot & tag == "medi"), 
            mapping = aes(x = x, y = y, fill = val)) +
  facet_wrap(~year, ncol = 1) +
  scale_fill_viridis_c(na.value = NA, "Prevalence") +
  # scale_fill_gradientn(colors = iddoPal::iddo_palettes$BlGyRd, 
  #                      "",
  #                      breaks = c(0, 0.5, 1), 
  #                      labels = c("0  (all K76)", "0.5", "1  (all 76T)"),
  #                      limits = c(0,1)) +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(title = "Median") +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.justification = "top")
p1

p2 <- ggplot() +
  geom_sf(data = afr, fill = "white") +
  geom_tile(data = df %>%
              filter(year %in% years_to_plot & tag == "sd"), 
            mapping = aes(x = x, y = y, fill = val)) +
  facet_wrap(~year, ncol = 1, strip.position = "right") +
  scale_fill_distiller(palette = "Oranges", 
                       na.value = NA, 
                       "Uncertainty", 
                       direction = 1,
                       trans = "sqrt") +
  xlab("Longitude") +
  ylab("") +
  labs(title = "Standard deviation") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.justification = "top")

pal <- iddoPal::iddo_palettes$soft_blues

df_sum <- df %>%
  filter(tag == "medi") %>%
  group_by(year) %>%
  summarise(q = list(quantile(val, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)))) %>%
  unnest_wider(q) %>%
  ungroup() %>%
  mutate(year = as.numeric(year))
# pivot_longer(cols = ends_with("%"),
#              names_to = "Quantile",
#              values_to = "val")



p3 <- ggplot(df_sum) +
  geom_line(aes(x = year, y = `0%`, linetype = "0% - 100%")) +
  geom_line(aes(x = year, y = `100%`, linetype = "0% - 100%")) +
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = "2.5% - 97.5%")) + #fill=pal[6]) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`, fill = "25% - 75%")) + #fill=pal[4]) +
  geom_ribbon(aes(x = year, ymin = `50%`, ymax = `50%`, fill = "50%")) + #fill=pal[1]) +
  geom_line(aes(x = year, y = `50%`), col = pal[1], linewidth = 1) +
  scale_linetype_manual("", values = c("0% - 100%" = 2)) +
  scale_fill_manual("", values = c("2.5% - 97.5%" = pal[6], "25% - 75%" = pal[4], "50%" = pal[1])) +
  ylab("Prevalence") +
  xlab("Year") +
  labs(title = "Estimated prevalence of Kelch 13 mutations in Africa") +
  theme_bw() +
  theme(legend.spacing.y = unit(-0.9, "cm"),
        legend.background = element_rect(fill = NA))
p3
ggsave("figures/k13_out_times.png", height = 2, width = 4, scale = 2)

# probably need to look at this next to data: 100% goes up before 2006
# snap box to extent of years/0

library(patchwork)

# would like to group x axis labels but giving up for now
# the feature should be available with axes/axis_titles arguments ....
# plot_space() is cool tho
p1 + plot_spacer() + p2 + plot_layout(ncol = 3, widths = c(4, -0.9, 4), guides = "collect")

ggsave("figures/k13_out.png", height = 3.6, width = 3, scale = 2.5)

layout <- "
AABBC#
AABBDD
"

# ecdf: surveillance PfPR?
q1 <- p1 + 
  p2 + 
  guide_area() + 
  plot_spacer() + 
  plot_layout(design = layout, guides = "collect")
q1

q1 + inset_element(p3, left = 0, bottom = 0, right = 1, top = 0.4)
library(gridExtra)






