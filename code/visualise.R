# do some visualisation in here ... 
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
    theme(legend.position = "bottom")
  
  p
}

survey_effort_panel(c("output/circmat_k13/surveillance_effort_k13.grd",
                      "output/circmat_crt/surveillance_effort_crt.grd",
                      "output/circmat_pfmdr86/surveillance_effort_pfmdr86.grd"),
                    lyr_names = c("Pfkelch13 (all)", "Pfcrt (76)", "Pfmdr1 (86)"))

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
        height = 6, width=15, scale=0.7)


# after all that, I might just facet_wrap ....



###############################################################################
# time !
pal <- iddoPal::iddo_palettes$soft_blues

pred_time_plot <- function(in_path, title=""){
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









