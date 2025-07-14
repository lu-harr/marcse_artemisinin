# do some visualisation in here ... 
setwd("~/Desktop/MARCSE/k13_seafrica")
source("code/setup.R")
library(viridisLite)
library(sf)

# mask
afr <- world %>%
  filter(continent == "Africa") %>%
  vect() %>%
  crop(ext(-21, 63, -35, 37)) %>%
  st_as_sf()

output_dir <- "output/circmat_pfmdr86/"

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

survey_effort_panel(c("output/circmat_k13/surveillance_effort_k13.grd",
                      "output/circmat_crt/surveillance_effort_crt.grd",
                      "output/circmat_pfmdr86/surveillance_effort_pfmdr86.grd"),
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


# after all that, I might just facet_wrap ....



###############################################################################
# time !

pred_time_plot <- function(in_path, 
                           title = "",
                           points_path = "",
                           pal = iddoPal::iddo_palettes$soft_blues){
  # wrapping up plot of all-Africa posterior meds over times
  preds <- rast(in_path)
  
  coords <- xyFromCell(preds, cells(preds))
  vals <- terra::extract(preds, coords)
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4),
           tag = substr(lyr, 11, 14)) %>% # pick out year and thingo
    filter(tag == "medi") %>%
    group_by(year) %>%
    summarise(q = list(quantile(val, c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)))) %>%
    unnest_wider(q) %>%
    ungroup() %>%
    mutate(year = as.numeric(year))
  
  p <- ggplot(df) +
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
    labs(title = title) +
    ylim(0, 1) +
    scale_x_continuous(breaks = seq(2000, 2022, 2), expand = c(0,0)) +
    theme_bw() +
    theme(legend.spacing.y = unit(-10, "cm"),
          legend.background = element_rect(fill = NA))
  
  if (points_path != ""){
    message("Watch out! I set limits manually!")
    mut_data <- setup_mut_data(points_path, min_year = 2000)
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


p1 <- pred_time_plot("output/circmat_k13/preds_all.grd",
               title = "(a) Pfkelch13")
p2 <- pred_time_plot("output/circmat_crt/preds_all.grd",
               title = "(b) Pfcrt K76T")
p3 <- pred_time_plot("output/circmat_pfmdr86/preds_all.grd",
               title = "(c) Pfmdr1 N86Y")
p4 <- pred_time_plot("output/circmat_pfmdr184/preds_all.grd",
               title = "(d) Pfmdr1 Y184F")
p5 <- pred_time_plot("output/circmat_pfmdr1246/preds_all.grd",
               title = "(e) Pfmdr1 D1246Y")

library(patchwork)
p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("figures/all_markers_time.png", scale = 1.5, height = 7, width = 6)

p1 <- pred_time_plot("output/circmat_k13/preds_all.grd",
                     title = "(a) Pfkelch13",
                     points_path = "data/moldm_k13_nomarker.csv")
p2 <- pred_time_plot("output/circmat_crt/preds_all.grd",
                     title = "(b) Pfcrt K76T",
                     points_path = "data/moldm_crt76.csv")
p3 <- pred_time_plot("output/circmat_pfmdr86/preds_all.grd",
                     title = "(c) Pfmdr1 N86Y",
                     points_path = "../moldm/clean/pfmdr_single_86.csv")
p4 <- pred_time_plot("output/circmat_pfmdr184/preds_all.grd",
                     title = "(d) Pfmdr1 Y184F",
                     points_path = "../moldm/clean/pfmdr_single_184.csv")
p5 <- pred_time_plot("output/circmat_pfmdr1246/preds_all.grd",
                     title = "(e) Pfmdr1 D1246Y",
                     points_path = "../moldm/clean/pfmdr_single_1246.csv")

p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 1, guides = "collect", axis_title = "collect")
ggsave("figures/all_markers_time_pts.png", scale = 1.5, height = 7, width = 6)


  

#########################################################################
# wrap up model predictions for pfcrt/pfmdr1:

library(RColorBrewer)
oranges <- brewer.pal(9, "Oranges")
blrd <- iddoPal::iddo_palettes$BlGyRd

map_pred_row <- function(in_path,
                         years,
                         pal,
                         field = c("medi", "sd"),
                         xlab = "Longitude",
                         ylab = "Latitude",
                         legend_lim = waiver(),
                         top_pan = FALSE){
  preds <- rast(in_path)
  coords <- xyFromCell(preds, cells(preds))
  vals <- terra::extract(preds, coords)
  df <- cbind(coords, vals) %>%
    pivot_longer(starts_with("2"),
                 names_to = "lyr",
                 values_to = "val") %>%
    mutate(year = substr(lyr, 1, 4),
           tag = substr(lyr, 11, 14)) %>%
    filter(year %in% years & tag == field) # pick out year and thingo
  
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
  
  if (field == "medi"){
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

years_to_plot <- c("2006","2010", "2014", "2018", "2022")
p1 <- map_pred_row("output/circmat_crt/preds_all.grd", 
             years = years_to_plot, field = "medi", pal = blrd,
             legend_lim = c(0,1), ylab = "", top_pan = TRUE)#
p2 <- map_pred_row("output/circmat_pfmdr86/preds_all.grd", 
                   years = years_to_plot, field = "medi", pal = blrd,
                   legend_lim = c(0,1), xlab = "", ylab = "", top_pan = TRUE)# ylab = "(a)")
p3 <- map_pred_row("output/circmat_pfmdr184/preds_all.grd", 
                   years = years_to_plot, field = "medi", pal = blrd,
                   legend_lim = c(0,1), xlab = "", ylab = "")#  ylab = "(c)")
p4 <- map_pred_row("output/circmat_pfmdr1246/preds_all.grd", 
                   years = years_to_plot, field = "medi", pal = blrd,
                   legend_lim = c(0,1), xlab = "",ylab = "")#  ylab = "(e)")
p5 <- map_pred_row("output/circmat_crt/preds_all.grd", 
                   years = years_to_plot, field = "sd", pal = oranges,
                   legend_lim = c(0,0.1), xlab = "",ylab = "")#  ylab = "(b)")
p6 <- map_pred_row("output/circmat_pfmdr86/preds_all.grd", 
                   years = years_to_plot, field = "sd", pal = oranges,
                   legend_lim = c(0,0.1), xlab = "",ylab = "")#  ylab = "(b)")
p7 <- map_pred_row("output/circmat_pfmdr184/preds_all.grd", 
                   years = years_to_plot, field = "sd", pal = oranges,
                   legend_lim = c(0,0.1), xlab = "",ylab = "")#  ylab = "(d)")
p8 <- map_pred_row("output/circmat_pfmdr1246/preds_all.grd", 
                   years = years_to_plot, field = "sd", pal = oranges,
                   legend_lim = c(0,0.1), xlab = "",ylab = "")#  ylab = "(f)")

library(gridExtra)
library(grid)
library(cowplot)


# chuck the rows together
pcol <- plot_grid(p2 + theme(legend.position = "none"), 
             p6 + theme(legend.position = "none"), 
             p3 + theme(legend.position = "none"),
             p7 + theme(legend.position = "none"),
             p4 + theme(legend.position = "none"),
             p8 + theme(legend.position = "none"),
             #ncol = 1, rel_heights = c(1.21, rep(1, 5))) +
             ncol = 1, rel_heights = c(1.167, rep(1, 5))) +
  theme(panel.spacing = unit(0, "cm"))

#pcol

#pcol <- grid.arrange(arrangeGrob(pcol, left = y.grob, bottom = x.grob))

legend1 = get_legend(p2)
legend2 = get_legend(p6)

plegend <- plot_grid(legend1, legend2, ncol = 1, axis = "r")

#plegend

p <- plot_grid(pcol, plegend, rel_widths = c(0.8,0.2))

# this df worked when I was using subfigure labels ...
# df <- data.frame(xmin = rep(0.02, 3),
#                  xmax = rep(0.062, 3),
#                  ymin = c(0.017, 0.34, 0.66),
#                  ymax = c(0.312, 0.635, 0.955),
#                  lab = c("D1246Y", "Y184F", "N86Y"))
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
                            label = lab), angle = 90) +
  geom_text(data = data.frame(x = rep(0.85, 2), y = c(0.34, 0.84), 
                              label = c("Estimate SD", "Prevalence")),
            aes(x = x, y = y, label = label))

p

# gave up on add_sub
ggsave("figures/mdr_out.png", height = 9, width = 8.5)

########################################################################
# can I fit crt in too?
p2 <- map_pred_row("output/circmat_pfmdr86/preds_all.grd", 
                   years = years_to_plot, field = "medi", pal = blrd,
                   legend_lim = c(0,1), xlab = "", ylab = "")# ylab = "(a)")

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

p

ggsave("figures/crt_mdr_out.png", height = 9, width = 7)







  
  