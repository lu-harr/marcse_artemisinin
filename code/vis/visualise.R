# do some visualisation in here ... 
setwd("~/Desktop/MARCSE/k13_seafrica")
source("code/setup.R")
source("code/vis/vis_funcs.R")
library(sf)

library(viridisLite)
library(RColorBrewer)
library(iddoPal)
oranges <- brewer.pal(9, "Oranges")
blrd <- iddoPal::iddo_palettes$BlGyRd
iddoblue <- iddo_palettes_discrete$iddo[1]
# taking some bright colours that don't coincide with viridis:
case_pal <- c("#E37210", iddoblue, "#c7047c")

# mask
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

#output_dir <- "output/circmat_pfmdr86/"

###############################################################################
# surveillance effort

survey_effort_panel <- function(in_path, 
                                agg_factor = 1, 
                                pan = "", 
                                main = "",
                                xlab = "Longitude",
                                ylab = "Latitude",
                                lyr_names = main){
  
  test_dens <- lapply(in_path, rast) %>%
    rast()
  if (agg_factor > 1){
    test_dens <- aggregate(test_dens, agg_factor)
  }
  
  # nasty formula to convert to tests per 100kmsq (ish)
  test_dens <- test_dens * 100 / (res(test_dens)[1] * 111) ** 2
  
  test_dens <- mask(test_dens, st_as_sf(afr))
  
  # surveil <- xyFromCell(test_dens, cell = cells(test_dens)) %>%
  #   as.data.frame() %>%
  #   cbind(unlist(extract(test_dens, cells(test_dens))))
  surveil <- cbind(xyFromCell(test_dens, cell = cells(test_dens)), 
                   extract(test_dens, cells(test_dens)))
  names(surveil) <- c("x", "y", lyr_names)
  
  surveil <- surveil %>%
    pivot_longer(lyr_names, names_to = "lyr", values_to = "effort") %>%
    mutate(lyr = factor(lyr, levels = lyr_names))
  
  p <- ggplot() +
    geom_sf(data = st_as_sf(afr), fill = "white") + # not showing anything in the background here ...
    geom_tile(data = surveil, aes(x, y, fill = effort)) +
    geom_sf(data = st_as_sf(afr), colour = "white", fill = NA) +
    # scale_fill_viridis_c(na.value = NA, bquote(atop("Tests per","~100"~km^2)), trans="sqrt") +
    scale_fill_viridis_c(na.value = NA, bquote("Tests per ~100"~km^2), trans="sqrt") +
    facet_wrap(~ lyr) +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(paste(pan, main)) +
    theme(legend.position = "bottom", strip.text.x = element_text(size = 12))
  
  p
}

survey_effort_panel(c("output/k13/surveillance_effort_k13.grd",
                      "output/crt76/surveillance_effort_crt76.grd",
                      "output/mdr86/surveillance_effort_mdr86.grd"),
                    lyr_names = c("Pfkelch13", "Pfcrt K76T", "Pfmdr1 N86Y"))

# p1 <- survey_effort_panel("output/circmat_k13/surveillance_effort_k13.grd", 
#                           pan = "(a)", main = "Kelch 13", xlab="")
# p2 <- survey_effort_panel("output/circmat_crt/surveillance_effort_crt.grd", 
#                           pan = "(b)", main = "Pfcrt 76", ylab="")
# p3 <- survey_effort_panel("output/pfmdr_hier/surveillance_effort_pfmdr86.grd", 
#                           pan = "(c)", main = "Pfmdr1 86/184/1246", xlab="", ylab="")
# 
# library(patchwork)
# p1 + p2 + p3 +
#   plot_layout(ncol = 3, widths = c(4, 4, 4),
#               axis_titles = "collect")
ggsave("figures/surveillance_effort_all.png", 
        height = 5, width=10, scale = 0.8)

survey_effort_panel(c("output/k13_marcse/surveillance_effort_k13.grd",
                      "output/crt76/surveillance_effort_crt76.grd",
                      "output/mdr86/surveillance_effort_mdr86.grd"),
                    lyr_names = c("Pfkelch13", "Pfcrt K76T", "Pfmdr1 N86Y"))
ggsave("figures/surveillance_effort_all_marcse.png", 
       height = 5, width=10, scale = 1.1)

survey_effort_panel(c("output/k13/surveillance_effort_k13.grd",
                      "output/k13_marcse/surveillance_effort_k13_marcse.grd"),
                    lyr_names = c("moldm", "moldm + unpublished data"))
ggsave("~/Desktop/presentations/MARCSE/surveillance_effort_marcse.png", 
       height = 5, width=8, scale = 0.9)

# after all that, I might just facet_wrap ....

###############################################################################
# time !
# note to self: this has the potential to crash ya laptop

p1 <- pred_time_plot("output/k13_marcse/gneiting_sparse/",
               title = "(a) Pfkelch13")
p2 <- pred_time_plot("output/crt76/gneiting_sparse/",
               title = "(b) Pfcrt K76T")
p3 <- pred_time_plot("output/mdr86/gneiting_sparse/",
               title = "(c) Pfmdr1 N86Y")
p4 <- pred_time_plot("output/mdr184/gneiting_sparse/",
               title = "(d) Pfmdr1 Y184F")
p5 <- pred_time_plot("output/mdr1246/gneiting_sparse/",
               title = "(e) Pfmdr1 D1246Y")

library(patchwork)
p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("figures/all_markers_time_bin.png", scale = 1.5, height = 7, width = 6)


p1 <- pred_time_plot("output/k13_marcse/bb_gne/",
                     title = "(a) Pfkelch13")
p2 <- pred_time_plot("output/crt76/bb_gne/",
                     title = "(b) Pfcrt K76T")
p3 <- pred_time_plot("output/mdr86/bb_gne/",
                     title = "(c) Pfmdr1 N86Y")
p4 <- pred_time_plot("output/mdr184/bb_gne/",
                     title = "(d) Pfmdr1 Y184F")
p5 <- pred_time_plot("output/mdr1246/bb_gne/",
                     title = "(e) Pfmdr1 D1246Y")

p1 + p2 + p3 + p4 + p5 + 
  plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("figures/all_markers_time_bb.png", scale = 1.5, height = 7, width = 6)

p1 <- pred_time_plot("output/k13_marcse/gneiting_ahmc/preds_all.tif",
                     title = "(a) Pfkelch13")

p1 + p1 + p1 + p1 + p1 + 
  plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("~/Desktop/presentations/test_time2.png", scale = 1.5, height = 7, width = 6)


# make a version of this with highlights to k13 panel, highlights to other pred panels?
# after all that, not sure this adds much?
zambezi <- list(xmin = 19, xmax = 26, ymin = -21, ymax = -14)
victoria <- list(xmin = 27, xmax = 35, ymin = -6, ymax = 2)
eswatini <- list(xmin = 34, xmax = 41, ymin = 12, ymax = 19)

zoom_df <- rbind(zambezi, victoria, eswatini) %>%
  as.data.frame() %>%
  unnest() %>%
  suppressWarnings()

# ended up abandoning this
# p1 <- pred_time_plot("output/k13/gneiting_sparse/preds_all.grd",
#                      title = "(a) Pfkelch13",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "data/clean/moldm_k13_nomarker.csv")
# p2 <- pred_time_plot("output/crt76/gneiting_sparse/preds_all.grd",
#                      title = "(b) Pfcrt K76T",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "data/clean/moldm_crt76.csv")
# p3 <- pred_time_plot("output/mdr86/gneiting_sparse/preds_all.grd",
#                      title = "(c) Pfmdr1 N86Y",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "../moldm/clean/pfmdr_single_86.csv")
# p4 <- pred_time_plot("output/mdr184/gneiting_sparse/preds_all.grd",
#                      title = "(d) Pfmdr1 Y184F",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "../moldm/clean/pfmdr_single_184.csv")
# p5 <- pred_time_plot("output/mdr1246/gneiting_sparse/preds_all.grd",
#                      title = "(e) Pfmdr1 D1246Y",
#                      zooms = zoom_df, zoom_pal = case_pal,
#                      points_path = "../moldm/clean/pfmdr_single_1246.csv")
# 
# p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
# ggsave("figures/all_markers_time_pts_zoom_gneiting.png", scale = 1.5, height = 7, width = 6)

#########################################################################
# wrap up model predictions for pfcrt/pfmdr1:

# make sure you've run what's up the top of this script
years_to_plot <- c("2004","2009", "2014", "2019", "2024")
p1 <- map_pred_row("output/crt76/bb_gne/preds_medians.tif", 
             years = years_to_plot, field = "50", pal = blrd,
             legend_lim = c(0,1), ylab = "", top_pan = TRUE)#
p2 <- map_pred_row("output/mdr86/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), xlab = "", ylab = "")# ylab = "(a)")
p3 <- map_pred_row("output/mdr184/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), xlab = "", ylab = "")#  ylab = "(c)")
p4 <- map_pred_row("output/mdr1246/bb_gne/preds_medians.tif", 
                   years = years_to_plot, field = "50", pal = blrd,
                   legend_lim = c(0,1), xlab = "",ylab = "")#  ylab = "(e)")
p5 <- map_pred_row("output/mdr86/bb_gne/preds_sds.tif", 
                   years = years_to_plot, field = "sd", pal = oranges,
                   legend_lim = c(0,0.2), xlab = "",ylab = "")#  ylab = "(b)")
p6 <- map_pred_row("output/mdr86/bb_gne/preds_sds.tif", 
                   years = years_to_plot, field = "sd", pal = oranges,
                   legend_lim = c(0,0.2), xlab = "",ylab = "")#  ylab = "(b)")
p7 <- map_pred_row("output/mdr184/bb_gne/preds_sds.tif", 
                   years = years_to_plot, field = "sd", pal = oranges,
                   legend_lim = c(0,0.2), xlab = "",ylab = "")#  ylab = "(d)")
p8 <- map_pred_row("output/mdr1246/bb_gne/preds_sds.tif", 
                   years = years_to_plot, field = "sd", pal = oranges,
                   legend_lim = c(0,0.2), xlab = "",ylab = "")#  ylab = "(f)")
# legend_lim was c(0, 0.2) for the main fig ... although this probably could have been less
library(gridExtra)
library(grid)
library(cowplot)


# # chuck the rows together
# pcol <- plot_grid(p2 + theme(legend.position = "none"), 
#              p6 + theme(legend.position = "none"), 
#              p3 + theme(legend.position = "none"),
#              p7 + theme(legend.position = "none"),
#              p4 + theme(legend.position = "none"),
#              p8 + theme(legend.position = "none"),
#              #ncol = 1, rel_heights = c(1.21, rep(1, 5))) +
#              ncol = 1, rel_heights = c(1.167, rep(1, 5))) +
#   theme(panel.spacing = unit(0, "cm"))
# 
# #pcol
# 
# #pcol <- grid.arrange(arrangeGrob(pcol, left = y.grob, bottom = x.grob))
# 
# legend1 = get_legend(p2)
# legend2 = get_legend(p6)
# 
# plegend <- plot_grid(legend1, legend2, ncol = 1, axis = "r")
# 
# #plegend
# 
# p <- plot_grid(pcol, plegend, rel_widths = c(0.8,0.2))
# 
# # this df worked when I was using subfigure labels ...
# # df <- data.frame(xmin = rep(0.02, 3),
# #                  xmax = rep(0.062, 3),
# #                  ymin = c(0.017, 0.34, 0.66),
# #                  ymax = c(0.312, 0.635, 0.955),
# #                  lab = c("D1246Y", "Y184F", "N86Y"))
# df <- data.frame(xmin = rep(0.001, 3),
#                  xmax = rep(0.028, 3),
#                  ymin = c(0.0155, 0.337, 0.663),
#                  ymax = c(0.315, 0.639, 0.965),
#                  lab = c("D1246Y", "Y184F", "N86Y"))
# p <- p + 
#   # tried adding outer margin but that didn't do anything
#   geom_rect(data = df, aes(xmin=xmin, xmax=xmax, ymin=ymin, 
#                           ymax=ymax), 
#             colour="grey10", fill="grey85", linewidth=0.3) +
#   geom_text(data = df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
#                             label = lab), angle = 90) +
#   geom_text(data = data.frame(x = rep(0.85, 2), y = c(0.34, 0.84), 
#                               label = c("Estimate SD", "Prevalence")),
#             aes(x = x, y = y, label = label))
# 
# p
# 
# # gave up on add_sub
# ggsave("figures/mdr_out_bb.png", height = 9, width = 8.5)

########################################################################
# can I fit crt in too?
# p2 <- map_pred_row("output/mdr86/bb_gne/preds_all.tif", 
#                    years = years_to_plot, field = "50", pal = blrd,
#                    legend_lim = c(0,1), xlab = "", ylab = "")# ylab = "(a)")

# chuck the rows together
pcol <- plot_grid(p1 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")), 
                  p5 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
                  p2 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")), 
                  p6 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")), 
                  p3 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
                  p7 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
                  p4 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
                  p8 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
                  #ncol = 1, rel_heights = c(1.21, rep(1, 5))) +
                  ncol = 1, rel_heights = c(1.22, rep(1, 7))) +
  theme(panel.spacing = unit(0, "cm"))

#pcol

#pcol <- grid.arrange(arrangeGrob(pcol, left = y.grob, bottom = x.grob))

legend1 = get_legend(p2)
legend2 = get_legend(p6)

plegend <- plot_grid(legend1, legend2, ncol = 1, axis = "r")

#plegend

p <- plot_grid(pcol, plegend, rel_widths = c(0.8,0.22))

# this df worked when I was using subfigure labels ...
# df <- data.frame(xmin = rep(0.02, 3),
#                  xmax = rep(0.062, 3),
#                  ymin = c(0.017, 0.34, 0.66),
#                  ymax = c(0.312, 0.635, 0.955),
#                  lab = c("D1246Y", "Y184F", "N86Y"))
df <- data.frame(xmin = rep(0.001, 4),
                 xmax = rep(0.035, 4),
                 ymin = c(0.005, 0.249, 0.492, 0.735),
                 ymax = c(0.242, 0.487, 0.73, 0.975),
                 lab = c("Pfmdr1 D1246Y", "Pfmdr1 Y184F", "Pfmdr1 N86Y", "Pfcrt K76T"))
p <- p + 
  # tried adding outer margin but that didn't do anything
  geom_rect(data = df, aes(xmin=xmin, xmax=xmax, ymin=ymin,
                           ymax=ymax),
            colour="grey10", fill="grey85", linewidth=0.3) +
  geom_text(data = df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2,
                           label = lab), angle = 90) +
  geom_text(data = data.frame(x = rep(0.85, 2), y = c(0.34, 0.84), 
                              label = c("Estimate SD", "Prevalence")),
            aes(x = x, y = y, label = label))

ggsave("figures/crt_mdr_out_bb.png", p, height = 9, width = 7)

##############################################################################
# need a short version for powerpoint

# chuck the rows together
# pcol <- plot_grid(p1 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")), 
#                   p2 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")), 
#                   p3 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   p4 + theme(legend.position = "none", plot.margin = unit(rep(0,4), "cm")),
#                   ncol = 1, rel_heights = c(1.22, rep(1, 3))) +
#   theme(panel.spacing = unit(0, "cm"))

# I'm exhausted by this
years_to_plot <- c(2006, 2014, 2022)
preds <- rast(paste0("output/", c("crt76", "mdr86", "mdr184", "mdr1246"), "/gneiting_sparse/preds_all.grd")) %>%
  aggregate(fact=3)
preds <- preds[[grep("median", names(preds))]]
names(preds) <- paste0(names(preds), "_", rep(c("crt76", "mdr86", "mdr184", "mdr1246"), each = 25))
coords <- xyFromCell(preds, cells(preds))
vals <- terra::extract(preds, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         tag = substr(lyr, 11, 14),
         marker = str_extract(lyr, "[^_]+$")) %>%
  filter(year %in% years_to_plot & tag == "medi") # pick out year and thingo

df <- df %>% mutate(marker = case_when(marker == "crt76" ~ "Pfcrt K76T",
                                       marker == "mdr86" ~ "Pfmdr1 N86Y",
                                       marker == "mdr184" ~ "Pfmdr1 Y184F",
                                       marker == "mdr1246" ~ "Pfmdr1 D1246Y")) %>%
  mutate(marker = factor(marker, levels = c("Pfcrt K76T", "Pfmdr1 N86Y", "Pfmdr1 Y184F", "Pfmdr1 D1246Y")))

ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") + 
  geom_tile(aes(x = x, y = y, fill = val), data = df) +
  scale_fill_gradientn(name = "Prevalence",
                       colors = blrd, 
                       breaks = c(0, 0.5, 1), 
                       labels = c("0  (all wildtype)", "0.5", "1  (all mutant)"),
                       limits = c(0,1)) +
  facet_grid(year ~ marker) +
  scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  xlab("Longitude") +
  ylab("Latitude") +
  theme_bw()

ggsave("~/Desktop/presentations/MARCSE/crt_mdr_out.png", scale = 1.7, height = 3, width = 5)



# I'm exhausted by this
years_to_plot <- c(2006, 2014, 2022)
preds <- rast(paste0("output/", c("crt76", "mdr86", "mdr184", "mdr1246"), "/gneiting_sparse/preds_all.grd")) %>%
  aggregate(fact=3)
preds <- preds[[grep("sd", names(preds))]]
names(preds) <- paste0(names(preds), "_", rep(c("crt76", "mdr86", "mdr184", "mdr1246"), each = 25))
coords <- xyFromCell(preds, cells(preds))
vals <- terra::extract(preds, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         tag = substr(lyr, 11, 14),
         marker = str_extract(lyr, "[^_]+$"))
df <- filter(df, year %in% years_to_plot) # pick out year and thingo

df <- df %>% mutate(marker = case_when(marker == "crt76" ~ "Pfcrt K76T",
                                       marker == "mdr86" ~ "Pfmdr1 N86Y",
                                       marker == "mdr184" ~ "Pfmdr1 Y184F",
                                       marker == "mdr1246" ~ "Pfmdr1 D1246Y")) %>%
  mutate(marker = factor(marker, levels = c("Pfcrt K76T", "Pfmdr1 N86Y", "Pfmdr1 Y184F", "Pfmdr1 D1246Y")))

ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") + 
  geom_tile(aes(x = x, y = y, fill = val), data = df) +
  scale_fill_gradientn(name = "Uncertainty",
                       colors = oranges, trans="sqrt") +
  facet_grid(year ~ marker) +
  scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  theme_bw() +
  xlab("Longitude") +
  ylab("Latitude")

ggsave("~/Desktop/presentations/MARCSE/crt_mdr_out_sd.png", scale = 1.7, height = 3, width = 5)

##############################################################################

coords <- xyFromCell(pfpr, cells(pfpr))
vals <- terra::extract(pfpr, coords)
df <- cbind(coords, vals) %>%
  mutate(cell = 1:nrow(coords)) %>%
  pivot_longer(starts_with("pfpr"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = as.numeric(substr(lyr, 6, 9))) # pick out year and thingo

cells_to_plot <- sample(1:nrow(coords), 50, replace = FALSE)
  
df_summ <- df %>%
  group_by(year) %>%
  summarise(q = list(quantile(val, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)))) %>%
  unnest_wider(q) %>%
  ungroup() %>%
  mutate(year = as.numeric(year))

pal = iddoPal::iddo_palettes$soft_blues

p <- ggplot(df_summ) +
  geom_line(aes(x = year, y = `0%`, linetype = "0% - 100%")) +
  geom_line(aes(x = year, y = `100%`, linetype = "0% - 100%")) +
  geom_ribbon(aes(x = year, ymin = `2.5%`, ymax = `97.5%`, fill = "2.5% - 97.5%")) + #fill=pal[6]) +
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`, fill = "25% - 75%")) + #fill=pal[4]) +
  geom_ribbon(aes(x = year, ymin = `50%`, ymax = `50%`, fill = "50%")) + #fill=pal[1]) +
  geom_line(aes(x = year, y = `50%`), col = pal[1], linewidth = 1) +
  scale_linetype_manual("", values = c("0% - 100%" = 2)) +
  scale_fill_manual("", values = c("2.5% - 97.5%" = pal[6], "25% - 75%" = pal[4], "50%" = pal[1])) +
  geom_line(data = df %>% filter(cell %in% cells_to_plot), 
            mapping = aes(x = year, y = val, colour = cell, group = cell)) +
  ylab("PfPR") +
  xlab("Year") +
  theme_bw() +
  theme(legend.spacing.y = unit(-10, "cm"),
        legend.background = element_rect(fill = NA))

p
  

##############################################################################
# compare obs to pred
# could also plot obs against prediction uncertainty

#mut_data$pred[is.na(mut_data$pred)] = 0.4

cex_transform <- function(from){
  sqrt(from) / sqrt(max(from))
  #log10(from) / log10(max(from)) *4
}




##############################################################################
# the same but ggplot

# gosh I prefer how base plotting deals with panels :/
# manipulate these further when I decide whether I want facetting etc
# obs_prev_panel("data/clean/moldm_k13_nomarker.csv",
#                "output/k13/circmat_sparse/preds_all.grd",
#                "k13 circmat", xlim = c(0, 0.4), ylim = c(0, 0.4))#,
#                #facet_bins = c(2010, 2014, 2021))

# obs_prev_panel("data/clean/moldm_k13_nomarker.csv",
#                 "output/k13/gneiting_sparse/preds_all.grd",
#                 "k13 gneiting", xlim = c(0, 0.4), ylim = c(0, 0.4))

obs_prev_panel("data/clean/moldm_marcse_k13_nomarker.csv",
               "output/k13_marcse/gneiting_ahmc/preds_all.tif",
               "", xlim = c(0, 0.6), ylim = c(0, 0.6), 
               ave_tag = "_50")
ggsave("figures/residuals_k13m_gne.png", height = 3, width = 4, scale = 2)

obs_prev_panel("data/clean/pfmdr_single_86.csv",
                "output/mdr86/gneiting_ahmc/preds_all.tif",
                "mdr86 gneiting", 
               ave_tag = "_50") #, facet_bins = c(2008, 2012, 2016, 2020))
ggsave("figures/residuals_86_gne.png", height = 9, width = 5, scale = 2)

obs_prev_panel("data/clean/pfmdr_single_86.csv",
               "output/mdr86/bb_gne/preds_all.tif",
               "mdr86 gneiting", 
               ave_tag = "_50")
ggsave("figures/residuals_86_bb.png", height = 9, width = 5, scale = 2)

# obs_prev_panel("data/clean/pfmdr_single_86.csv",
#                "output/mdr86/circmat/preds_all.grd", "mdr86 circmat", 
#                ave_tag = "_post_median") 
#                #facet_bins = c(2008, 2012, 2016, 2020) )
# ggsave("~/Desktop/residuals_86.png", height = 9, width = 5, scale = 2)

obs_prev_panel("data/clean/pfmdr_single_1246.csv",
                "output/mdr1246/gneiting_ahmc/preds_all.tif",
               #facet_bins = c(2008, 2012, 2016, 2020),
                "mdr1246 gneiting", 
               ave_tag = "_50")
ggsave("figures/residuals_1246_gne.png", height = 9, width = 5, scale = 2)

obs_prev_panel("data/clean/pfmdr_single_1246.csv",
               "output/mdr1246/bb_gne/preds_all.tif",
               #facet_bins = c(2008, 2012, 2016, 2020),
               "mdr1246 gneiting", buffer = 100000)
ggsave("figures/residuals_1246_bb.png", height = 9, width = 5, scale = 2)

# obs_prev_panel("data/clean/pfmdr_single_1246.csv",
#                "output/mdr1246/circmat/preds_all.grd",
#                #facet_bins = c(2008, 2012, 2016, 2020),
#                "mdr1246 circmat", 
#                ave_tag = "_post_median")
# ggsave("~/Desktop/residuals_1246.png", height = 9, width = 5, scale = 2)

obs_prev_panel("data/clean/pfmdr_single_184.csv",
                "output/mdr184/gneiting_ahmc/preds_all.tif",
                "mdr184 gneiting", 
               ave_tag = "_50")
ggsave("figures/residuals_184_gne.png", height = 9, width = 5, scale = 2)

# oops got rid of this
# obs_prev_panel("data/clean/pfmdr_single_184.csv",
#                "output/mdr184/circmat/preds_all.grd",
#                 "mdr184 circmat", 
#                ave_tag = "_post_median")

obs_prev_panel("data/clean/pfmdr_single_184.csv",
               "output/mdr184/bb_gne/preds_all.tif",
               "mdr184 gneiting", buffer = 100000)
ggsave("figures/residuals_184_bb.png", height = 9, width = 5, scale = 2)

# how do I capture clustered data where the model picks somewhere in the middle
# but there's lots of variance at the same location?
# rolling window?
# this is where the CV would be helpful to make my point ...

##############################################################################
# visualise coverages

coverages_fig("output/k13_marcse/gneiting_sparse/")
coverages_fig("output/k13_marcse/bb_gne/")
coverages_fig("output/crt76/gneiting_sparse/")
coverages_fig("output/crt76/bb_gne/")
coverages_fig("output/mdr184/gneiting_sparse/")
coverages_fig("output/mdr184/bb_gne/")
#coverages_fig("output/mdr86/gneiting_sparse/")
#coverages_fig("output/mdr86/bb_gne/")
coverages_fig("output/mdr1246/gneiting_sparse/")
coverages_fig("output/mdr1246/bb_gne/")

# here's what I put in the MS
p1 <- coverages_fig(list("output/k13_marcse/gneiting_sparse/", 
                         "output/k13_marcse/bb_gne/"))
p2 <- coverages_fig(list("output/mdr184/gneiting_sparse/", 
                         "output/mdr184/bb_gne/"))

p <- plot_grid(p1 + theme(legend.position = "none"), 
               p2 + theme(legend.position = "none"), nrow = 2, 
               labels = "auto", label_x = 0.1, label_y = 0.97)
plot_grid(p, get_legend(p1), rel_widths = c(1, 0.25))
ggsave("figures/coverages.png", height = 5, width = 7)

mut_data <- read_rds("output/k13_marcse/gneiting_sparse/mut_data.rds")
lower <- rast("output/k13_marcse/gneiting_sparse/preds_lower.tif")
upper <- rast("output/k13_marcse/gneiting_sparse/preds_upper.tif")

lower_upper_panel <- function(path, 
                           main = "", 
                           show_nas = FALSE, 
                           pal = colorRamp(iddo_palettes$BlGyRd),
                           xlim = c(0,1), # define limits to pred/obs panel
                           ylim = c(0,1), # define limits to pred/obs panel
                           facet_bins = NULL, # apply facets over time?
                           ave_tag = "_50", # mean? median? what are the surfaces called in the stack?
                           buffer = 1, # option to reland points?
                           bb = NULL){
  
  mut_data <- read_rds(paste0(path, "mut_data.rds")) %>%
    arrange(present/tested)
  lower <- rast(paste0(path, "preds_lower.tif"))
  upper <- rast(paste0(path, "preds_upper.tif"))
  yrs_pred <- str_extract(names(upper), "\\d{4}")
  
  # get predictions for each row in `mut_data`
  mut_data$lower <- NA
  mut_data$upper <- NA
  yrs_to_extract <- unique(mut_data$year)
  for (yr in yrs_to_extract){
    if (yr %in% yrs_pred){
      idx <- which(mut_data$year == yr)
      vall <- terra::extract(lower[[paste0(yr, "_2.5")]], 
                            mut_data[idx, c("x", "y")],
                            ID = FALSE, search_radius = buffer)
      valu <- terra::extract(upper[[paste0(yr, "_97.5")]], 
                            mut_data[idx, c("x", "y")],
                            ID = FALSE, search_radius = buffer)
      if(ncol(vall) < 3){
        # idk why we have to have an inconsistent return when |idx| == 1
        mut_data[idx, "lower"] <- vall[1,1]
        mut_data[idx, "upper"] <- valu[1,1]
      } else{
        mut_data[idx, "lower"] <- vall[, paste0(yr, "_2.5")]
        mut_data[idx, "upper"] <- valu[, paste0(yr, "_97.5")]
      }
    }
  }
  
  message(nrow(mut_data))
  mut_data <- mut_data %>% 
    filter(!is.na(lower) & !is.na(upper)) %>%
    mutate(idx = 1:nrow(.), 
           covered = ifelse(present/tested > lower & present/tested < upper, TRUE, FALSE))
  message(nrow(mut_data))
  
  mut_data_long <- mut_data %>%
    pivot_longer(cols = c(lower, upper), names_to = "lu", values_to = "value")
  
  ggplot(mut_data) +
    geom_point(aes(x = idx, y = present/tested,
                   size = tested, 
                   col = covered), 
               pch = 1) +
    geom_line(data = mut_data_long, aes(x = idx, y = value, group = idx, col = covered)) +
    scale_colour_manual(values = iddo_palettes$iddo, "CI covers\nobservation") +
    scale_size_continuous(range = c(0.1, 6), "Tested") +
    xlab("Index") +
    ylab("Observed prevalence") +
    theme_bw() #+
    # for some reason, this legend.position doesn't play with get_legend
    # theme(legend.position = "bottom")
    
}


lower_upper_panel("output/k13_marcse/gneiting_sparse/")
lower_upper_panel("output/k13_marcse/bb_gne/")
lower_upper_panel("output/mdr86/gneiting_sparse/")
lower_upper_panel("output/mdr86/bb_gne/")
# not 100% sure what's going on with crt ... may need refitting ..
lower_upper_panel("output/crt76/gneiting_sparse/")
lower_upper_panel("output/crt76/bb_gne/")
lower_upper_panel("output/mdr1246/gneiting_sparse/")
lower_upper_panel("output/mdr1246/bb_gne/")
lower_upper_panel("output/mdr184/gneiting_sparse/")
lower_upper_panel("output/mdr184/bb_gne/")


p1 <- lower_upper_panel("output/k13_marcse/gneiting_sparse/")
p2 <- lower_upper_panel("output/k13_marcse/bb_gne/")
p3 <- lower_upper_panel("output/mdr184/gneiting_sparse/")
p4 <- lower_upper_panel("output/mdr184/bb_gne/")

legs <- get_legend(p1)
# get_plot_component(p1, "guide-box", return_all = TRUE)

p <- plot_grid(p1 + theme(legend.position = "none", axis.title = element_blank()), 
          p2 + theme(legend.position = "none", axis.title = element_blank()), 
          p3 + theme(legend.position = "none", axis.title = element_blank()), 
          p4 + theme(legend.position = "none", axis.title = element_blank()),
          nrow = 4, align = "v", labels = "auto", label_y = 0.95, label_x = 0.05)
plot_grid(p, legs, nrow = 1, rel_widths = c(1, 0.1))
ggsave("figures/crints.png", scale = 1, height = 7, width = 9)
# check this again with zeroes removed
# and have a look at upper and lower bound surfaces?

##############################################################################
# a plot of all preds in all years

preds <- rast("output/k13_marcse/gneiting_sparse/preds_all.grd")
coords <- xyFromCell(preds, cells(preds))
vals <- terra::extract(preds, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         tag = substr(lyr, 11, 14),
         marker = str_extract(lyr, "[^_]+$")) %>%
  filter(tag == "medi") # pick out year and thingo

ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") + 
  geom_tile(aes(x = x, y = y, fill = val), data = df) +
  scale_fill_viridis_c(na.value = NA, "Prevalence", trans = "sqrt") +
  facet_wrap(~year, ncol = 6) +
  scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  xlab("Longitude") +
  ylab("Latitude")

ggsave("~/Desktop/presentations/MARCSE/k13_all_years.png", scale = 1.7, height = 4, width = 5)

preds <- rast("output/k13_marcse/bb_gne/preds_all.tif")
coords <- xyFromCell(preds, cells(preds))
vals <- terra::extract(preds, coords)
df <- cbind(coords, vals) %>%
  pivot_longer(starts_with("2"),
               names_to = "lyr",
               values_to = "val") %>%
  mutate(year = substr(lyr, 1, 4),
         tag = substr(lyr, 11, 14),
         marker = str_extract(lyr, "[^_]+$")) %>%
  filter(marker == "50") # pick out year and thingo

ggplot() +
  geom_sf(data = st_as_sf(afr), fill = "white") + 
  geom_tile(aes(x = x, y = y, fill = val), data = df) +
  scale_fill_viridis_c(na.value = NA, "Prevalence", trans = "sqrt") +
  facet_wrap(~year, ncol = 6) +
  scale_x_continuous(breaks = seq(-20, 40, 20), "Longitude") +
  scale_y_continuous(breaks = seq(-20, 40, 20), "Latitude") +
  xlab("Longitude") +
  ylab("Latitude")

ggsave("~/Desktop/presentations/MARCSE/k13_all_years_bb.png", scale = 1.7, height = 4, width = 5)



